#!/usr/bin/env python3
"""Build maize-pathogen classification benchmark: 573 positives + 500 curated negatives."""

import json
import os
import re
import time
from collections import defaultdict

import requests
from Bio import SeqIO

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
OUT = os.path.join(BASE, "docs", "validation", "classification_benchmark")
POS_FASTA = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")
CATALOG = json.load(open("/Users/gfgao/Desktop/blacksoil_metaG/Figshare/taxonomy.json"))
NEG_FASTA = os.path.join(OUT, "negative_queries.fasta")
NEG_META = os.path.join(OUT, "negative_queries_meta.json")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"
TOOL = "maize_pathogen_db"
DELAY = 0.4

CAT_TARGETS = {"bacteria": 100, "viruses": 100, "fungi": 250, "oomycetes": 50}
MARKER = {
    "bacteria": '("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title] OR "16S ribosomal RNA gene"[Title])',
    "viruses": '("complete genome"[Title] OR "genome"[Title] OR "polyprotein"[Gene])',
    "fungi": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title] OR "18S ribosomal RNA"[Gene])',
    "oomycetes": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title] OR "18S ribosomal RNA"[Gene])',
}

NEGATIVE_TAXA = {
    "bacteria": [
        "Streptomyces coelicolor", "Streptomyces lividans", "Streptomyces griseus",
        "Streptomyces avermitilis", "Streptomyces venezuelae", "Bacillus subtilis",
        "Bacillus amyloliquefaciens", "Pseudomonas putida", "Rhizobium leguminosarum",
        "Bradyrhizobium japonicum", "Sinorhizobium meliloti", "Nitrosomonas europaea",
        "Nitrobacter winogradskyi", "Rhodococcus erythropolis", "Arthrobacter globiformis",
        "Micrococcus luteus", "Sphingomonas paucimobilis", "Caulobacter crescentus",
        "Methylobacterium extorquens", "Rhodobacter sphaeroides", "Paracoccus denitrificans",
        "Geobacter sulfurreducens", "Shewanella oneidensis", "Deinococcus radiodurans",
        "Thermus thermophilus", "Lactobacillus plantarum", "Lactococcus lactis",
        "Bifidobacterium longum", "Corynebacterium glutamicum", "Mycobacterium smegmatis",
        "Lactobacillus acidophilus", "Lactobacillus rhamnosus", "Lactobacillus casei",
        "Lactobacillus reuteri", "Lactobacillus helveticus", "Lactobacillus delbrueckii",
        "Leuconostoc mesenteroides", "Oenococcus oeni", "Pediococcus acidilactici",
        "Streptococcus thermophilus", "Clostridium acetobutylicum", "Clostridium beijerinckii",
        "Clostridium cellulolyticum", "Clostridium thermocellum", "Ruminococcus flavefaciens",
        "Ruminococcus albus", "Fibrobacter succinogenes", "Prevotella ruminicola",
        "Bacteroides thetaiotaomicron", "Akkermansia muciniphila", "Faecalibacterium prausnitzii",
        "Roseburia intestinalis", "Eubacterium rectale", "Bifidobacterium adolescentis",
        "Bifidobacterium bifidum", "Bifidobacterium breve", "Bifidobacterium animalis",
        "Collinsella aerofaciens", "Zymomonas mobilis", "Gluconobacter oxydans",
        "Acetobacter aceti", "Acetobacter pasteurianus", "Komagataeibacter xylinus",
        "Methylococcus capsulatus", "Methylosinus trichosporium", "Nitrosospira multiformis",
        "Nitrosococcus oceani", "Nitrospira defluvii", "Nitrospira moscoviensis",
        "Desulfovibrio vulgaris", "Desulfovibrio desulfuricans", "Desulfobacter postgatei",
        "Desulfotomaculum reducens", "Syntrophobacter fumaroxidans", "Pelobacter carbinolicus",
        "Geobacter metallireducens", "Geobacter uraniireducens", "Anaeromyxobacter dehalogenans",
        "Myxococcus xanthus", "Stigmatella aurantiaca", "Sorangium cellulosum",
        "Chondromyces crocatus", "Herpetosiphon aurantiacus", "Chloroflexus aurantiacus",
        "Roseiflexus castenholzii", "Dehalococcoides mccartyi", "Synechocystis sp. PCC 6803",
        "Synechococcus elongatus", "Anabaena sp. PCC 7120", "Nostoc punctiforme",
        "Microcystis aeruginosa", "Prochlorococcus marinus", "Trichodesmium erythraeum",
        "Lactobacillus fermentum", "Lactobacillus salivarius", "Lactobacillus paracasei",
        "Lactobacillus johnsonii", "Lactobacillus gasseri", "Lactobacillus crispatus",
        "Lactobacillus jensenii", "Propionibacterium freudenreichii", "Amycolatopsis mediterranei",
        "Saccharopolyspora erythraea", "Salinispora arenicola", "Streptomyces clavuligerus",
        "Streptomyces hygroscopicus", "Streptomyces albus", "Frankia alni",
    ],
    "fungi": [
        "Mortierella elongata", "Mortierella alpina", "Chaetomium globosum",
        "Humicola grisea", "Humicola fuscoatra", "Podospora anserina",
        "Neurospora crassa", "Neurospora sitophila", "Aspergillus nidulans",
        "Aspergillus oryzae", "Penicillium rubens", "Penicillium citrinum",
        "Penicillium roqueforti", "Trichoderma reesei", "Trichoderma longibrachiatum",
        "Acremonium chrysogenum", "Stachybotrys chartarum", "Metarhizium anisopliae",
        "Beauveria bassiana", "Cordyceps militaris", "Clonostachys rosea",
        "Thielavia terrestris", "Myceliophthora thermophila", "Thermomyces lanuginosus",
        "Talaromyces marneffei", "Aspergillus terreus", "Penicillium camemberti",
        "Penicillium commune", "Penicillium paneum", "Mucor racemosus",
        "Lichtheimia corymbifera", "Cunninghamella elegans", "Umbelopsis ramanniana",
        "Komagataella phaffii", "Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
        "Yarrowia lipolytica", "Pichia kudriavzevii", "Debaryomyces hansenii",
        "Kluyveromyces lactis", "Blastobotrys adeninivorans", "Scheffersomyces stipitis",
        "Sugiyamaella lignohabitans", "Wickerhamomyces anomalus", "Rhodotorula toruloides",
        "Cutaneotrichosporon oleaginosum", "Apiotrichum porosum", "Lentinula edodes",
        "Pleurotus ostreatus", "Agaricus bisporus", "Coprinopsis cinerea",
        "Schizophyllum commune", "Trametes versicolor", "Ganoderma lucidum",
        "Phanerochaete chrysosporium", "Trichoderma atroviride", "Chaetomium thermophilum",
        "Thermothelomyces thermophilus", "Sordaria macrospora", "Podospora curvicolla",
        "Mortierella parvispora", "Mortierella hyalina", "Mortierella verticillata",
        "Mucor hiemalis", "Backusella circina", "Absidia glauca", "Actinomucor elegans",
        "Pirella circinans", "Radiomyces spectabilis", "Syncephalastrum racemosum",
        "Thamnidium elegans", "Circinella muscae", "Gongronella butleri",
        "Cunninghamella bertholletiae", "Metarhizium robertsii", "Isaria fumosorosea",
        "Purpureocillium lilacinum", "Paecilomyces variotii", "Lecanicillium lecanii",
        "Simplicillium lanosoniveum", "Tolypocladium inflatum", "Engyodontium album",
        "Scopulariopsis brevicaulis", "Microascus brevicaulis", "Cephalotheca foveolata",
        "Coniochaeta ligniaria", "Diatrype disciformis", "Xylaria hypoxylon",
        "Daldinia concentrica", "Hypoxylon fragiforme", "Biscogniauxia nummularia",
        "Nemania serpens", "Coprinellus disseminatus", "Coprinus comatus",
        "Agaricus subrufescens", "Volvariella volvacea", "Pleurotus eryngii",
        "Flammulina velutipes", "Hypsizygus marmoreus", "Pholiota nameko",
        "Grifola frondosa", "Ganoderma sinense", "Trametes hirsuta", "Trametes pubescens",
        "Fomes fomentarius", "Fomitopsis pinicola", "Postia placenta", "Gloeophyllum trabeum",
        "Irpex lacteus", "Phlebia radiata", "Ceriporiopsis subvermispora",
        "Dichomitus squalens", "Bjerkandera adusta", "Pycnoporus cinnabarinus",
        "Stereum hirsutum", "Auricularia auricula-judae", "Tremella fuciformis",
        "Meyerozyma guilliermondii", "Ogataea polymorpha", "Cyberlindnera jadinii",
        "Kazachstania exigua", "Torulaspora delbrueckii", "Zygosaccharomyces bailii",
        "Zygosaccharomyces rouxii", "Hanseniaspora uvarum", "Metschnikowia pulcherrima",
        "Lachancea thermotolerans", "Saccharomycodes ludwigii", "Starmerella bombicola",
        "Moesziomyces antarcticus", "Sporidiobolus salmonicolor", "Rhodotorula glutinis",
        "Rhodosporidium toruloides", "Cystobasidium minutum", "Leucosporidium creatinivorum",
    ],
    "oomycetes": [
        "Saprolegnia parasitica", "Saprolegnia diclina", "Achlya bisexualis",
        "Achlya klebsiana", "Leptolegnia chapmanii", "Thraustotheca clavata",
        "Aphanomyces astaci", "Aphanomyces invadans", "Halophytophthora vesicula",
        "Apodachlya pyrifera", "Leptomitus lacteus", "Dictyuchus sterilis",
        "Saprolegnia ferax", "Saprolegnia hypogyna", "Achlya americana", "Achlya conspicua",
        "Achlya racemosa", "Achlya flagellata", "Achlya colorata", "Dictyuchus monosporus",
        "Dictyuchus pseudodictyum", "Isoachlya anisospora", "Isoachlya toruloides",
        "Leptolegnia caudata", "Aphanodictyon papillatum", "Plectospira myriandra",
        "Calyptralegnia achlyoides", "Brevilegnia unisperma", "Brevilegnia megasperma",
        "Sapromyces elongatus", "Apodachlya brachynema", "Halophytophthora polymorphica",
        "Halophytophthora batemanensis", "Halophytophthora avicennae", "Halophytophthora exoprolifera",
    ],
    "viruses": [
        "Escherichia phage Lambda", "Enterobacteria phage T4", "Enterobacteria phage T7",
        "Bacillus phage phi29", "Mycobacterium phage D29", "Pseudomonas phage phi6",
        "Synechococcus phage S-PM2", "Prochlorococcus phage P-SSM2", "Vibrio phage VHML",
        "Lactococcus phage c2", "Staphylococcus phage 2638A", "Salmonella phage P22",
        "Acinetobacter phage phiAB1", "Burkholderia phage phiE125", "Caulobacter phage phiCbK",
        "Autographa californica nucleopolyhedrovirus", "Bombyx mori nucleopolyhedrovirus",
        "Cotesia congregata bracovirus", "Amsacta moorei entomopoxvirus", "Choristoneura fumiferana NPV",
        "Enterobacteria phage T1", "Enterobacteria phage T5", "Enterobacteria phage phiX174",
        "Enterobacteria phage MS2", "Enterobacteria phage Qbeta", "Enterobacteria phage PRD1",
        "Enterobacteria phage PM2", "Enterobacteria phage M13", "Enterobacteria phage f1",
        "Enterobacteria phage fd", "Enterobacteria phage P2", "Enterobacteria phage P4",
        "Enterobacteria phage Mu", "Enterobacteria phage N4", "Enterobacteria phage 186",
        "Enterobacteria phage HK97", "Enterobacteria phage HK022", "Escherichia phage P1",
        "Enterobacteria phage RB69", "Enterobacteria phage RB49", "Bacillus phage SPP1",
        "Bacillus phage SPO1", "Bacillus phage PBS1", "Streptococcus phage Cp-1",
        "Streptococcus phage Dp-1", "Lactococcus phage 936", "Lactococcus phage P335",
        "Listeria phage A511", "Listeria phage P100", "Pseudomonas phage phiKZ",
        "Pseudomonas phage PAK_P1", "Pseudomonas phage PP7", "Pseudomonas phage phi8",
        "Salmonella phage epsilon15", "Salmonella phage FelixO1", "Vibrio phage ICP1",
        "Vibrio phage CTXphi", "Synechococcus phage P60", "Synechococcus phage S-CRM01",
        "Prochlorococcus phage MED4-213", "Prochlorococcus phage P-SSM4", "Cyanophage Ma-LMM01",
        "Mycobacterium phage TM4", "Mycobacterium phage Bxb1", "Mycobacterium phage Che9c",
        "Mycobacterium phage L5", "Helicoverpa zea nucleopolyhedrovirus",
        "Spodoptera frugiperda nucleopolyhedrovirus", "Plutella xylostella granulovirus",
        "Cydia pomonella granulovirus", "Adoxophyes orana nucleopolyhedrovirus",
        "Epiphyas postvittana nucleopolyhedrovirus", "Penicillium chrysogenum virus",
        "Saccharomyces cerevisiae virus L-A",
    ],
}


def normalize_seq(seq):
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


def binomials(name):
    out = set()
    for part in re.split(r"[/,()]", name):
        m = re.match(r"\s*([A-Za-z]+)\s+([a-z][a-z-]+)", part)
        if m:
            out.add((m.group(1).lower(), m.group(2).lower()))
    return out


def resolve_taxid(name):
    r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi",
                     params={"db": "taxonomy", "term": f'"{name}"[Organism]', "retmode": "json",
                             "email": EMAIL, "tool": TOOL}, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None, None
    tid = ids[0]
    r2 = requests.get(f"{ENTREZ_BASE}/esummary.fcgi",
                      params={"db": "taxonomy", "id": tid, "retmode": "json",
                              "email": EMAIL, "tool": TOOL}, timeout=30)
    r2.raise_for_status()
    data = r2.json()["result"]
    uid = data["uids"][0]
    return tid, data[uid].get("scientificname", "")


def fetch_sequences(taxid, cat, retmax=10):
    query = f'txid{taxid}[Organism] AND {MARKER[cat]}'
    r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi",
                     params={"db": "nucleotide", "term": query, "retmax": retmax, "retmode": "json",
                             "email": EMAIL, "tool": TOOL}, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    r2 = requests.get(f"{ENTREZ_BASE}/efetch.fcgi",
                      params={"db": "nucleotide", "id": ",".join(ids), "rettype": "fasta", "retmode": "text",
                              "email": EMAIL, "tool": TOOL}, timeout=60)
    r2.raise_for_status()
    out = []
    for block in r2.text.strip().split("\n\n"):
        lines = [ln.strip() for ln in block.split("\n") if ln.strip()]
        if len(lines) < 2 or not lines[0].startswith(">"):
            continue
        seq = "".join(lines[1:])
        out.append({"header": lines[0][1:], "seq": seq})
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    pos_seqs = set()
    for rec in SeqIO.parse(POS_FASTA, "fasta"):
        pos_seqs.add(normalize_seq(str(rec.seq)))

    catalog_binomials = set()
    catalog_taxids = set()
    for r in CATALOG:
        catalog_taxids.add(str(r["taxid"]))
        catalog_binomials |= binomials(r["species"])

    negative = []
    done_taxids = set()
    if os.path.exists(NEG_META) and os.path.exists(NEG_FASTA):
        existing_meta = json.load(open(NEG_META))
        seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(NEG_FASTA, "fasta")}
        negative = [{**m, "seq": seqs[m["qid"]]} for m in existing_meta if m["qid"] in seqs]
        done_taxids = {m["taxid"] for m in negative}
        print(f"Resuming with {len(negative)} existing negatives", flush=True)
    for cat, names in NEGATIVE_TAXA.items():
        print(f"Downloading {cat} negatives...", flush=True)
        for name in names:
            if len([x for x in negative if x["category"] == cat]) >= CAT_TARGETS[cat]:
                break
            try:
                taxid, resolved = resolve_taxid(name)
                time.sleep(DELAY)
                if not taxid or taxid in catalog_taxids or taxid in done_taxids:
                    continue
                done_taxids.add(taxid)
                for rec in fetch_sequences(taxid, cat):
                    time.sleep(DELAY)
                    norm = normalize_seq(rec["seq"])
                    if len(norm) < 100 or norm in pos_seqs:
                        continue
                    cand_bins = binomials(resolved)
                    if cand_bins & catalog_binomials:
                        continue
                    qid = f"n{len(negative) + 1:04d}|{taxid}|{cat}"
                    negative.append({
                        "qid": qid, "taxid": taxid, "species": resolved,
                        "category": cat, "accession": rec["header"].split()[0],
                        "title": rec["header"], "seq": rec["seq"],
                    })
                    if len([x for x in negative if x["category"] == cat]) >= CAT_TARGETS[cat]:
                        break
            except Exception as e:
                print(f"  {name}: {e}", flush=True)
            print(f"  {name} -> {sum(1 for x in negative if x['category'] == cat)}/{CAT_TARGETS[cat]}", flush=True)

    with open(NEG_FASTA, "w") as f:
        for rec in negative:
            f.write(f">{rec['qid']} {rec['title']}\n")
            for i in range(0, len(rec["seq"]), 80):
                f.write(rec["seq"][i:i + 80] + "\n")
    meta = [{k: rec[k] for k in ("qid", "taxid", "species", "category", "accession", "title")}
            for rec in negative]
    with open(NEG_META, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    from collections import Counter
    print("\nDownloaded negatives:")
    for cat, target in CAT_TARGETS.items():
        n = sum(1 for x in negative if x["category"] == cat)
        print(f"  {cat}: {n}/{target}")
    print("Total:", len(negative))


if __name__ == "__main__":
    main()
