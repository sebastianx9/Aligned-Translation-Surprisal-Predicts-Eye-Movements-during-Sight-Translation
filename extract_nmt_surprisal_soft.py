"""
Compute both argmax-aligned and soft-aligned NMT surprisal.

Soft-aligned surprisal for source word w:
  c_soft(w) = sum_t  s_t * ( word_attn[t,w] / sum_j word_attn[t,j] )

i.e. each Czech token's surprisal is distributed across all source words
proportionally to (normalised) cross-attention weight.
"""

import csv, math, os, torch
import numpy as np
from transformers import MarianMTModel, MarianTokenizer

SENTENCES_CSV  = "/Users/sebastianx/eyetracked-multi-modal-translation/probes/Sentences.csv"
OUTPUT_CSV     = "/Users/sebastianx/Dissertation_Data/nmt_surprisal_soft_word.csv"
MODEL_NAME     = "Helsinki-NLP/opus-mt-en-cs"
DEVICE         = torch.device("cpu")

print(f"Loading {MODEL_NAME} …")
tokenizer = MarianTokenizer.from_pretrained(MODEL_NAME)
model     = MarianMTModel.from_pretrained(MODEL_NAME, output_attentions=True)
model.eval(); model.to(DEVICE)

def load_sentences(path):
    sentences = {}
    with open(path, newline="", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            sid, text = line.split(",", 1)
            sentences[sid.strip()] = text.strip()
    return sentences

def subwords_to_word_map(tokens):
    word_map, wi = [], -1
    for tok in tokens:
        if tok in ("<pad>", "</s>", "<unk>"):
            word_map.append(-1); continue
        if tok.startswith("▁") or wi == -1:
            wi += 1
        word_map.append(wi)
    return word_map

def process_sentence(text):
    enc = tokenizer([text], return_tensors="pt", truncation=True, max_length=128).to(DEVICE)
    src_tokens = tokenizer.convert_ids_to_tokens(enc["input_ids"][0].tolist())

    with torch.no_grad():
        gen_ids = model.generate(enc["input_ids"],
                                 attention_mask=enc["attention_mask"],
                                 num_beams=4, max_length=200)[0]

    decoder_input_ids = gen_ids[:-1].unsqueeze(0)
    target_ids        = gen_ids[1:].unsqueeze(0)

    with torch.no_grad():
        out = model(input_ids=enc["input_ids"],
                    attention_mask=enc["attention_mask"],
                    decoder_input_ids=decoder_input_ids,
                    output_attentions=True)

    log_probs = torch.nn.functional.log_softmax(out.logits[0], dim=-1)
    surprisals = [-log_probs[t, tok].item()
                  for t, tok in enumerate(target_ids[0].tolist())]

    attn_mean = out.cross_attentions[-1][0].mean(dim=0).cpu().numpy()  # (T-1, S)

    cs_tokens  = tokenizer.convert_ids_to_tokens(target_ids[0].tolist())
    en_words   = text.split()
    en_word_map = subwords_to_word_map(src_tokens)
    n_en = len(en_words)
    n_cs = len(cs_tokens)

    # word-level attention matrix (n_cs, n_en)
    word_attn = np.zeros((n_cs, n_en))
    for sub_idx, wi in enumerate(en_word_map):
        if 0 <= wi < n_en and sub_idx < attn_mean.shape[1]:
            word_attn[:, wi] += attn_mean[:, sub_idx]

    # argmax and soft surprisal accumulators
    hard_surp = {i: [] for i in range(n_en)}
    soft_surp = np.zeros(n_en)

    for t, (tok, surp) in enumerate(zip(cs_tokens, surprisals)):
        if tok in ("</s>", "<pad>", "<unk>"): continue
        row = word_attn[t]
        row_sum = row.sum()

        # argmax
        hard_surp[int(row.argmax())].append(surp)

        # soft (normalised distribution)
        if row_sum > 0:
            soft_surp += surp * (row / row_sum)

    rows = []
    for wi, word in enumerate(en_words):
        hs = hard_surp[wi]
        rows.append({
            "word_index":        wi,
            "word":              word,
            "surprisal_sum":     round(sum(hs), 6) if hs else None,   # argmax
            "surprisal_soft":    round(float(soft_surp[wi]), 6),       # soft
        })
    return rows

sentences = load_sentences(SENTENCES_CSV)
print(f"Loaded {len(sentences)} sentences.\n")

all_rows = []
for sid in sorted(sentences.keys()):
    text = sentences[sid]
    print(f"  {sid}: {text[:60]}")
    try:
        rows = process_sentence(text)
        for r in rows:
            r["sentence_id"] = sid
        all_rows.extend(rows)
    except Exception as e:
        print(f"    ERROR: {e}")

fields = ["sentence_id", "word_index", "word", "surprisal_sum", "surprisal_soft"]
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(all_rows)

print(f"\nDone. {len(all_rows)} rows → {OUTPUT_CSV}")
