
import argparse

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from transformers import AutoModelForCausalLM, AutoTokenizer
from tqdm.auto import tqdm


class TextDataset(Dataset):
    def __init__(self, filename, tokenizer, max_length=128):
        with open(filename, "r", encoding="utf-8") as f:
            texts = [line.strip() for line in f if line.strip()]

        if not texts:
            raise ValueError("Training file contains no text.")

        self.samples = tokenizer(
            texts,
            truncation=True,
            max_length=max_length,
            padding="max_length",
            return_tensors="pt",
        )

    def __len__(self):
        return len(self.samples["input_ids"])

    def __getitem__(self, index):
        return {
            "input_ids": self.samples["input_ids"][index],
            "attention_mask": self.samples["attention_mask"][index],
        }


def compute_kd_loss(
    teacher_logits,
    student_logits,
    input_ids,
    attention_mask,
    alpha=0.5,
    temperature=2.0,
):
    # Shift logits and labels for next-token prediction.
    s_logits = student_logits[:, :-1, :]
    t_logits = teacher_logits[:, :-1, :]
    labels = input_ids[:, 1:]

    # Ignore padded target tokens.
    mask = attention_mask[:, 1:].bool()

    # Cross-entropy loss
    ce = F.cross_entropy(
        s_logits.reshape(-1, s_logits.size(-1)),
        labels.reshape(-1),
        reduction="none",
    )
    ce = ce.reshape_as(labels)
    ce_loss = ce[mask].mean()

    # Knowledge distillation loss: KL(teacher || student)
    T = temperature

    student_log_probs = F.log_softmax(
        s_logits.float() / T, dim=-1
    )
    teacher_log_probs = F.log_softmax(
        t_logits.float() / T, dim=-1
    )

    kd = F.kl_div(
        student_log_probs,
        teacher_log_probs,
        reduction="none",
        log_target=True,
    ).sum(dim=-1)

    kd_loss = kd[mask].mean() * (T ** 2)

    total_loss = alpha * ce_loss + (1 - alpha) * kd_loss

    return total_loss, ce_loss, kd_loss


def train(args):
    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    print("Device:", device)

    # Load tokenizer (same token IDs for teacher and student).
    tokenizer = AutoTokenizer.from_pretrained(args.teacher)

    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # Load models.
    teacher = AutoModelForCausalLM.from_pretrained(
        args.teacher
    ).to(device)

    student = AutoModelForCausalLM.from_pretrained(
        args.student
    ).to(device)

    # Full-logit KD requires matching vocabulary/token IDs.
    if teacher.config.vocab_size != student.config.vocab_size:
        raise ValueError("Teacher and student vocabularies differ.")

    teacher.eval()
    teacher.requires_grad_(False)

    student.train()

    dataset = TextDataset(
        args.data,
        tokenizer,
        max_length=args.max_length,
    )

    loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
    )

    optimizer = torch.optim.AdamW(
        student.parameters(),
        lr=args.learning_rate,
    )

    # Training loop
    for epoch in range(args.epochs):
        total = 0.0

        progress = tqdm(
            loader,
            desc=f"Epoch {epoch + 1}/{args.epochs}",
        )

        for batch in progress:
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)

            optimizer.zero_grad()

            # Teacher: forward only, no gradients.
            with torch.no_grad():
                t_outputs = teacher(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                )

            # Student: gradients enabled.
            s_outputs = student(
                input_ids=input_ids,
                attention_mask=attention_mask,
            )

            loss, ce_loss, kd_loss = compute_kd_loss(
                teacher_logits=t_outputs.logits,
                student_logits=s_outputs.logits,
                input_ids=input_ids,
                attention_mask=attention_mask,
                alpha=args.alpha,
                temperature=args.temperature,
            )

            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                student.parameters(), max_norm=1.0
            )
            optimizer.step()

            total += loss.item()

            progress.set_postfix(
                loss=f"{loss.item():.4f}",
                ce=f"{ce_loss.item():.4f}",
                kd=f"{kd_loss.item():.4f}",
            )

        avg_loss = total / len(loader)
        print(f"Epoch {epoch + 1}: avg loss = {avg_loss:.4f}")

    # Save fine-tuned student.
    student.save_pretrained(args.output)
    tokenizer.save_pretrained(args.output)

    print("Student saved to:", args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--teacher",
        default="Qwen/Qwen2.5-70B",
    )
    parser.add_argument(
        "--student",
        default="Qwen/Qwen2.5-32B",
    )
    parser.add_argument("--data", default="train.txt")
    parser.add_argument("--output", default="./student_kd")
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--batch_size", type=int, default=2)
    parser.add_argument("--max_length", type=int, default=128)
    parser.add_argument("--learning_rate", type=float, default=2e-5)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--temperature", type=float, default=2.0)

    args = parser.parse_args()
    train(args)
