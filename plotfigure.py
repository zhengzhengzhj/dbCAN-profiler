import pandas as pd
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context("poster")
def arg_parse():
    parser = argparse.ArgumentParser(description='plotting figure')
    parser.add_argument('-ab',"--abundance",help='abundance.',type=str)
    parser.add_argument('fig', help='which figures to plot.')
    parser.add_argument('-i', "--input",help='which figures to plot.')
    parser.add_argument('-t', "--type",default="class",help='which figures to plot.')
    parser.add_argument('-s1', "--sample1",help='which figures to plot.')
    parser.add_argument('-s2', "--sample2",help='which figures to plot.')
    parser.add_argument('-s1l', "--s1label",default="P1",help='which figures to plot.')
    parser.add_argument('-s2l', "--s2label",default="Golden")
    parser.add_argument('-l', "--level",default="class")
    return parser.parse_args()

def PlotTaxonomy(args):
    df = pd.read_csv(args.input,sep="\t")
    name = args.type
    #print (df[name])
    sns.countplot(x = name, data = df)
    sns.despine(ax=None, top=True, right=True, left=False, bottom=False)
    plt.xticks(rotation=20,fontsize=5)
    plt.show()

def barplot_Taxonomy(args):
    dfs = []
    for filename in open(args.input):
        filename = filename.rstrip("\n")
        df = pd.read_csv(filename,sep="\t")
        df["samples"] = filename.split(".")[0]
        dfs.append(df)
    dftables = pd.concat(dfs)
    #sns.countplot(x="samples",data=dftables,hue="class")
    sns.countplot(x="samples",data=dftables,hue=args.level)
    plt.show()


def Venn_fam(args):
    from matplotlib_venn import venn2, venn2_circles
    sample1 = pd.read_csv(args.sample1)  
    sample2 = pd.read_csv(args.sample2)  
    CAzy1 = set(sample1["CAZy"])
    CAzy2 = set(sample2["CAZy"])
    venn2([CAzy1,CAzy2],set_labels=(args.s1label,args.s2label))
    plt.show()

def filter_out_enzyme_number(table):
    bools = []
    for i in table["CAZy"]:
        if i[0].isdigit():
            bools.append(False)
        elif i in ["PL0","GH0","GT0","CBM0","AA0","CE0"]:
            bools.append(False)
        else:
            bools.append(True)
    table = table[bools]
    #print (table)
    return table

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#newstring = string.translate(table)
import re
def add_column_5type(table):
    cols = []
    for i in table["CAZy"]:
        fam = re.sub(r'[0-9]+', '', i)
        cols.append(fam)
        #print (i,fam)
    table["fam"] = cols
    return table

def FPKM_diff(args):
    sample1 = pd.read_csv(args.sample1)
    sample1["sample"] = "sample1"
    sample2 = pd.read_csv(args.sample2)  
    sample2["sample"] = "sample2"
    sample1 = filter_out_enzyme_number(sample1)
    sample2 = filter_out_enzyme_number(sample2)
    samples = pd.merge(sample1,sample2,on=["CAZy"],how="outer")
    samples.fillna(0,inplace=True)
    samples["diff_abs"] = np.abs(samples["FPKM_x"] - samples["FPKM_y"])
    samples["diff"] = samples["FPKM_x"] - samples["FPKM_y"]
    samples.rename(columns={"FPKM_x": args.s1label, "FPKM_y": args.s2label},inplace=True)
    samples.sort_values("diff_abs",inplace=True)

    # ### adding columns
    # add_column_5type(table)

    ax = samples.plot(x="CAZy", y=[args.s1label,args.s2label], kind="bar",figsize=(9,8))
    plt.xticks(rotation=90)
    axins = inset_axes(ax,  "30%", "30%" ,loc="upper center", borderpad=1)
    mean_value = np.mean(samples["diff_abs"])
    

    sns.kdeplot(samples["diff_abs"],shade=True,ax=axins,label=f"Mean:{round(mean_value,2)}")
    

    leg = plt.legend(frameon=False,handletextpad=-2.0, handlelength=0)
    for item in leg.legendHandles:
        item.set_visible(False)
    plt.show()



def servenal_kle(table,axins):

    uniqc = table["fam"].unique()
    for fam in uniqc:
        mean = round(np.mean(table[table["fam"] == fam]["diff_abs"]),2)
        sns.kdeplot(data=table[table["fam"] == fam],x="diff_abs" ,shade=True,ax=axins,label=fam + ": " + str(mean) ) 
    plt.legend(fontsize=14,frameon=False)
    #print (uniqc)


def classify_x(locs,labels):
    tmplabel =[re.sub(r'[0-9]+', '', i.get_text()) for i in labels]
    labelspos = {}
    for i in range(len(tmplabel)):
        labelspos.setdefault(tmplabel[i],[]).append(locs[i])

    tmppos = {}   
    for lab in labelspos:
        tmppos[lab] = (min(labelspos[lab]),max(labelspos[lab]))
    #print (tmppos)
    #return tmppos
    locations = [];labels = []
    for lab in tmppos:
        minx = tmppos[lab][0]
        maxx = tmppos[lab][1]
        labels.append("")
        locations.append(tmppos[lab][0])
        ##
        labels.append(lab)
        locations.append(minx/2+maxx/2)
        ##
        labels.append("")
        locations.append(tmppos[lab][1])
    return labels,locations

def set_xtick():
    locs, labels = plt.xticks()
    #print (locs)
    #print (labels)
    #plt.xticks([100],["MKMK"],fontsize=20)
    labels,locations = classify_x(locs,labels)
    plt.xticks(locations,labels,fontsize=20)
    #for i in labels:
    #    print (str(i))
    #plt.show()

colors_map =  {"GT":"darkorange","GH":"red","AA":"chartreuse","PL":"blue","CBM":"cyan"}


def set_xcolor(ax):
    asaaa = ax.get_xticklabels()
    for a in asaaa:
        fam = a.get_text()
        if fam in colors_map:
            cor = colors_map[fam]
            a.set_color(cor)

def Compare_diff(args):
    sample1 = pd.read_csv(args.sample1)
    sample1["sample"] = "sample1"
    sample2 = pd.read_csv(args.sample2)  
    sample2["sample"] = "sample2"
    sample1 = filter_out_enzyme_number(sample1)
    sample2 = filter_out_enzyme_number(sample2)
    samples = pd.merge(sample1,sample2,on=["CAZy"],how="outer")
    samples.fillna(0,inplace=True)
    samples["diff_abs"] = np.abs(samples["FPKM_x"] - samples["FPKM_y"])
    samples["diff"] = samples["FPKM_x"] - samples["FPKM_y"]
    samples.rename(columns={"FPKM_x": args.s1label, "FPKM_y": args.s2label},inplace=True)
    samples.sort_values(["CAZy","diff_abs"],inplace=True)
    #print (samples)
    ### adding columns
    samples = add_column_5type(samples)
    #print (samples)


    ### main figure
    ax = samples.plot(x="CAZy", y=[args.s1label,args.s2label], kind="bar",figsize=(10,8))
    plt.legend(frameon=False)
    plt.xticks(rotation=0)
    set_xtick()
    set_xcolor(ax)

    #plt.axis("off")
    axins = inset_axes(ax,  "30%", "30%" ,loc="upper center", borderpad=1)
  
    
    
    servenal_kle(samples,axins)
    #sns.kdeplot(data=samples,x="diff_abs" ,hue="fam", shade=True,ax=axins)    
    plt.xlim(-10,max(samples["diff_abs"]))
    #leg = plt.legend(frameon=False)
    #for item in leg.legendHandles:
    #    item.set_visible(False)
    

    ### veen 
    from matplotlib_venn import venn2, venn2_circles
    axins2 = inset_axes(ax,  "30%", "30%" ,loc="center left", borderpad=1)
    #sample1 = pd.read_csv(args.sample1)  
    #sample2 = pd.read_csv(args.sample2)  
    CAzy1 = set(sample1["CAZy"])
    CAzy2 = set(sample2["CAZy"])
    venn2([CAzy1,CAzy2],set_labels=(args.s1label,args.s2label),ax=axins2)

    plt.show()

def Fscore(args):
    data = pd.read_csv(args.input)
    abundance = pd.read_csv(args.abundance)
    data = data.merge(abundance,on=["otuid"])
    lineage_level = data[args.level].values
    #print (lineage_level)
    uniq_lineages = set(lineage_level)
    #print (data)
    df = data.groupby(args.level).sum()

    #df = pd.DataFrame()
    #for uniq_lineage in uniq_lineages:
    TP = df["TP"] ; TN = df["TN"] ; FP = df["FP"] ; FN = df["FN"]
    df["Acc"] = (df["TP"]+df["TN"])/(df["TP"]+df["TN"]+df["FP"]+df["FN"])
    df["Sen"] = TP/(TP+FN)
    df["Pre"] = TP/(TP+FP)
    df["Rec"] = TP/(TP+FN)
    df["F1"]  = 2*df["Pre"]*df["Rec"]/(df["Pre"]+df["Rec"])
    df["TPR"] = TP/(TP+FN)
    df["FPR"] = FP/(FP+TN)
    df[args.level] = df.index.values
    fig, ax = plt.subplots()
    for col in ["Acc","Sen","Pre","Rec","F1","relative_abund"]:
        df[col].plot(style="-o")
    plt.legend()
    plt.show()
    #print (data["otuid"])



def main():
    args = arg_parse()
    if args.fig == "Taxonomy":
        # python3 plotfigure.py Taxonomy -i None.Taxonomy
        PlotTaxonomy(args)
    if args.fig == "barplot_Taxonomy":
        # python3 plotfigure.py barplot_Taxonomy -i Taxonomy.filename
        barplot_Taxonomy(args)
    if args.fig == "Venn_fam":
        # python3 plotfigure.py Venn_fam -s1 -s2
        Venn_fam(args)
    if args.fig == "FPKM_diff":
        # python3 plotfigure.py FPKM_diff -s1 -s2
        FPKM_diff(args)
    if args.fig == "diff":
        # python3 plotfigure.py diff -s1 -s2
        Compare_diff(args)
    if args.fig == "Fscore":
        # python3 plotfigure.py Fscore -i out.csv -l c -abund sample1.abund
        Fscore(args)
if __name__ == "__main__":
    main()
