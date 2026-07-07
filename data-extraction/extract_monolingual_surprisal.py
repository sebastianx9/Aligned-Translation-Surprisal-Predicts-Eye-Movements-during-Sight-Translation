"""
Compute word-level monolingual surprisal for English source sentences
using GPT-2 (causal language model).

surprisal(word_t) = sum of -log P(subword_i | all previous tokens)
                    for all subwords in word_t

First word's surprisal = -log P from unconditional GPT-2 distribution
(no prior context), treated as valid but noted in output (n_tokens still set).

Output columns:
  sentence_id, word_index, word,
  surprisal_sum,   # sum of -log P for subwords in this word
  surprisal_mean,  # mean of -log P for subwords in this word
  n_tokens         # number of GPT-2 subword tokens in this word
"""

import argparse
import csv
import torch
from transformers import GPT2LMHeadModel, GPT2TokenizerFast


def load_sentences(path):
    sentences = {}
    with open(path, newline="", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            sid, text = line.split(",", 1)
            sentences[sid.strip()] = text.strip()
    return sentences


def compute_surprisal(sentence_text, tokenizer, model):
    enc = tokenizer(
        sentence_text,
        return_tensors="pt",
        return_offsets_mapping=True,
    )
    input_ids      = enc["input_ids"]
    offset_mapping = enc["offset_mapping"][0]

    with torch.no_grad():
        logits = model(input_ids).logits[0]

    log_probs = torch.nn.functional.log_softmax(logits, dim=-1)

    with torch.no_grad():
        logits_empty = model(
            torch.tensor([[tokenizer.eos_token_id]])
        ).logits[0, 0]
    lp0 = torch.nn.functional.log_softmax(logits_empty, dim=-1)

    ids = input_ids[0].tolist()
    token_surprisals = [-lp0[ids[0]].item()]
    for t in range(1, len(ids)):
        token_surprisals.append(-log_probs[t - 1, ids[t]].item())

    words = sentence_text.split()
    word_starts = []
    pos = 0
    for w in words:
        word_starts.append(pos)
        pos += len(w) + 1

    word_surprisals = {i: [] for i in range(len(words))}
    for tok_idx, (start, end) in enumerate(offset_mapping.tolist()):
        if start == end:
            continue
        surp = token_surprisals[tok_idx]
        wi = 0
        for k in range(len(word_starts) - 1, -1, -1):
            if start >= word_starts[k]:
                wi = k
                break
        word_surprisals[wi].append(surp)

    rows = []
    for wi, word in enumerate(words):
        sl = word_surprisals[wi]
        rows.append({
            "word_index":     wi,
            "word":           word,
            "surprisal_sum":  round(sum(sl), 6)           if sl else None,
            "surprisal_mean": round(sum(sl) / len(sl), 6) if sl else None,
            "n_tokens":       len(sl),
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description="Extract GPT-2 monolingual surprisal.")
    parser.add_argument("--sentences", required=True,
                        help="Path to Sentences.csv from the EMMT corpus")
    parser.add_argument("--output", required=True,
                        help="Output CSV path")
    args = parser.parse_args()

    print("Loading GPT-2 …")
    tokenizer = GPT2TokenizerFast.from_pretrained("gpt2")
    model     = GPT2LMHeadModel.from_pretrained("gpt2")
    model.eval()
    print("  Model loaded.\n")

    sentences = load_sentences(args.sentences)
    print(f"Loaded {len(sentences)} sentences.\n")

    all_rows = []
    for sid in sorted(sentences.keys()):
        text = sentences[sid]
        print(f"  {sid}: {text[:65]}")
        try:
            rows = compute_surprisal(text, tokenizer, model)
            for r in rows:
                r["sentence_id"] = sid
            all_rows.extend(rows)
        except Exception as e:
            print(f"    ERROR: {e}")

    fields = ["sentence_id", "word_index", "word",
              "surprisal_sum", "surprisal_mean", "n_tokens"]
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nDone. {len(all_rows)} word rows → {args.output}")


if __name__ == "__main__":
    main()
