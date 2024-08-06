import pymongo
import pandas as pd
client = pymongo.MongoClient('mongodb://127.0.0.1:30000')
db = client.GeneAge
db.authenticate('GeneAge_admin','GeneAge_admin_access',mechanism='SCRAM-SHA-1')




##conserveµÄ»ùÒò
with open("/home/mwshi/all_cluater", "r") as f:
    header = f.readline()
    header = header.replace('"','')
    w = open("/home/mwshi/conserve_age_result.csv", "w")
    w.write("Human gene age,Human gene age range,Mouse gene age,Mouse gene age range," + header)
    header = header.strip().split(",")
    for line in f:
        line = line.replace('"','')
        fields = line.strip().split(",")
        record = dict(zip(header, fields))
        human_age = db.Homo_sapiens.find_one({"ensembl_gene_id":record["Gene.stable.ID"]})
        mouse_age = db.Mus_musculus.find_one({"ensembl_gene_id":record["Mouse.gene.stable.ID"]})
        if human_age and mouse_age:
            w.write(",".join([str(human_age["gene_age"]), human_age["gene_Interval"],
                              str(mouse_age["gene_age"]), mouse_age["gene_Interval"]
                              ])+","+line)
        else:
            print(record["Gene.stable.ID"] + " and " + record["Mouse.gene.stable.ID"] + " not found")
    w.close()

