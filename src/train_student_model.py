
import argparse
import json
import torch

from datasets import Dataset
from transformers import AutoModelForCausalLM, AutoTokenizer
from peft import LoraConfig
from trl import SFTTrainer, SFTConfig


def load_training_data(path):
    records = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            sample = json.loads(line)

            if sample.get("verified") is not True:
                continue

            records.append({
                "prompt": [
                    {
                        "role": "user",
                        "content": sample["prompt"],
                    }
                ],
                "completion": [
                    {
                        "role": "assistant",
                        "content": sample["completion"],
                    }
                ],
            })

    if not records:
        raise ValueError("No verified training examples.")

    return Dataset.from_list(records)


def main(args):
    dataset = load_training_data(args.data)

    tokenizer = AutoTokenizer.from_pretrained(
        args.model
    )

    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    dtype = (
        torch.bfloat16
        if torch.cuda.is_available()
        and torch.cuda.is_bf16_supported()
        else torch.float32
    )

    model = AutoModelForCausalLM.from_pretrained(
        args.model,
        torch_dtype=dtype,
    )

    lora_config = LoraConfig(
        r=16,
        lora_alpha=32,
        lora_dropout=0.05,
        bias="none",
        task_type="CAUSAL_LM",
        target_modules=[
            "q_proj", "k_proj", "v_proj", "o_proj",
            "gate_proj", "up_proj", "down_proj"
        ],
    )

    training_config = SFTConfig(
        output_dir=args.output,
        num_train_epochs=args.epochs,
        per_device_train_batch_size=1,
        gradient_accumulation_steps=4,
        learning_rate=2e-4,
        max_length=args.max_length,
        completion_only_loss=True,
        gradient_checkpointing=True,
        bf16=dtype == torch.bfloat16,
        logging_steps=10,
        save_strategy="epoch",
        report_to="none",
    )

    trainer = SFTTrainer(
        model=model,
        args=training_config,
        train_dataset=dataset,
        processing_class=tokenizer,
        peft_config=lora_config,
    )

    trainer.train()

    trainer.save_model(args.output)
    tokenizer.save_pretrained(args.output)

    print("Trained LoRA adapter saved:", args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--model",
        default="Qwen/Qwen2.5-Coder-0.5B-Instruct"
    )
    parser.add_argument(
        "--data", default="data/verified.jsonl"
    )
    parser.add_argument(
        "--output", default="./qwen_kd"
    )
    parser.add_argument(
        "--epochs", type=int, default=3
    )
    parser.add_argument(
        "--max_length", type=int, default=2048
    )

    main(parser.parse_args())
