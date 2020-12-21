import csv
import argparse

parser = argparse.ArgumentParser(description = 'Script to classify infection', 
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input', type = str, dest = 'input', help = 'Input file containing mapping table', required = True)
parser.add_argument('--output', type = str, dest = 'output', help = 'Output file path', required = True)
parser.add_argument('--bias', type = float, dest = 'bias', help = 'Output file path', default = 5000.0, required = False)
parser.add_argument('--targets', type = int, dest = 'targets', help = 'Output file path', default = 3, required = False)
parser.add_argument('--target_names', type = str, nargs = "+", dest = 'target_names', help = 'Output file path', 
                                      default = ("ORF_1a", "5p_spike_cds", 
                                                 "O-linked_glycan_residue_region", "charite-similar_gene_E"), 
                                      required = False)
parser.add_argument('--sample-sheet', type = str, dest = 'samples', help = 'Sample file used to identify pools', required = True)


args, unknown = parser.parse_known_args()

def read_sample_sheet(path):
    samples = dict()
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter = ",")
        for line in reader:
            sample_id = line["sample_name"]
            if(sample_id in samples):
                raise Exception("Sample ID used twice: %s"%(sample_id,))
            samples[sample_id] = (line["i7"], line["i5"])
    return samples

def predict(data, bias_factor, targets, pred_function):
    corr_data = dict()
    for pool in data:
        corr_data[pool] = dict()
        for t in targets:
            bias = sum([int(data[pool][s][t]) for s in data[pool]]) / bias_factor
            for s in data[pool]:
                if(s not in corr_data[pool]):
                    corr_data[pool][s] = dict()
                corr_data[pool][s][t] = int(data[pool][s][t]) - bias
    prediction = dict()
    for pool in corr_data:
        prediction[pool] = dict()
        for s in corr_data[pool]:
            prediction[pool][s] = pred_function(corr_data[pool][s])
    return prediction

samples = read_sample_sheet(args.samples)

data = dict()

with open(args.input, "r") as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for line in reader:
        sample_id = line["id"]
        pool = samples[sample_id]
        if(pool not in data):
            data[pool] = dict()
        data[pool][sample_id] = line

result = predict(data, args.bias, args.target_names, lambda x: (sum([x[k] > 0 for k in x]) >= args.targets))

with open(args.output, "w") as f:
    writer = csv.DictWriter(f, fieldnames=["id", "status"], delimiter = "\t")
    writer.writeheader()
    for pool in result:
        for sample in result[pool]:
            writer.writerow({"id": sample, "status": result[pool][sample]})


