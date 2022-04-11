#import all pacakges

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import anndata

import seaborn as sns

from scipy.stats import rankdata 
from scipy.stats.stats import pearsonr   

import matplotlib.patches as ptc
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.axis import Tick 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#set path to object
date="20220321"
path_to_fast="/rds/general/user/aalmeid2"
path_out=path_to_fast+"/projects/cardiac_single_cell_biology/live/DCM_project/"
path_in=path_to_fast+"/projects/cardiac_single_cell_biology/live/DCM_project/"

adata=sc.read_h5ad(path_in + "ec_20210804_no_Doublets.h5ad")

def gene_coloc(adata,Gene1="INHBA",Gene2="NRG1",groupby=None,threshold=0.1,mode=("plot","table"),color1="Reds",color2="Greens",color3="Blues",plot_only_coloc=False,alpha=1,size=0.3,fontsize=8,label_width="25%",label_height="2.5%"):
    import warnings
    warnings.filterwarnings("ignore")
    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    warnings.filterwarnings("ignore")   
    def plot_coloc():
        df=pd.DataFrame(adata.obsm["X_umap"],columns=["x","y"],index=adata.obs.index)
        df["x1"]=df["x"]
        df["y1"]=df["y"]
        df["x2"]=df["x"]
        df["y2"]=df["y"]
        df["x3"]=df["x"]
        df["y3"]=df["y"]
        df["x4"]=df["x"]
        df["y4"]=df["y"]
        matrix["Gene1"]=matrix[Gene1]
        matrix["Gene2"]=matrix[Gene2]
        coloc=(matrix["Gene1"]+matrix["Gene2"])/2
        df.x1[matrix.Gene1<threshold]="NA"
        df.y1[matrix.Gene1<threshold]="NA"
        df.x1[matrix.Gene2>=threshold]="NA"
        df.y1[matrix.Gene2>=threshold]="NA"
        df.x2[matrix.Gene2<threshold]="NA"
        df.y2[matrix.Gene2<threshold]="NA"
        df.x2[matrix.Gene1>=threshold]="NA"
        df.y2[matrix.Gene1>=threshold]="NA"
        df.x3[matrix.Gene1>=threshold]="NA"
        df.y3[matrix.Gene1>=threshold]="NA"
        df.x3[matrix.Gene2>=threshold]="NA"
        df.y3[matrix.Gene2>=threshold]="NA"
        df.x4[matrix.Gene1<threshold]="NA"
        df.y4[matrix.Gene1<threshold]="NA"
        df.x4[matrix.Gene2<threshold]="NA"
        df.y4[matrix.Gene2<threshold]="NA"
        df.x1=df.x1.replace("NA",np.NaN)
        df.y1=df.y1.replace("NA",np.NaN)
        df.x2=df.x2.replace("NA",np.NaN)
        df.y2=df.y2.replace("NA",np.NaN)
        df.x3=df.x3.replace("NA",np.NaN)
        df.y3=df.y3.replace("NA",np.NaN)
        df.x4=df.x4.replace("NA",np.NaN)
        df.y4=df.y4.replace("NA",np.NaN)
        per_exp=df.x4.count()*100/(df.x1.count()+df.x2.count()+df.x4.count())
        per_total=df.x4.count()*100/(df.x1.count()+df.x2.count()+df.x4.count()+df.x3.count())
        per_gene1=df.x4.count()*100/(df.x1.count()+df.x4.count())
        per_gene2=df.x4.count()*100/(df.x2.count()+df.x4.count())
        c1Big = cm.get_cmap(color1, 512)
        c1 = ListedColormap(c1Big(np.linspace(0.3, 1.0, 256)))
        c2Big = cm.get_cmap(color2, 512)
        c2 = ListedColormap(c2Big(np.linspace(0.3, 1.0, 256)))
        c3Big = cm.get_cmap(color3, 512)
        c3 = ListedColormap(c3Big(np.linspace(0.3, 1.0, 256)))
        fig, ax = plt.subplots()
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False,labelleft=False)
        ax.scatter(df.x3,df.y3,s=size,c="lightgrey",alpha=1,marker="o",edgecolors="none")
        ax.set_title(Gene1 +"/" + Gene2 + " Co-localisation")
        if plot_only_coloc==True:
            im=ax.scatter(df.x1,df.y1,s=size,c="lightgrey",alpha=1,marker="o",edgecolors="none")
            im2=ax.scatter(df.x2,df.y2,s=size,c="lightgrey",alpha=1,marker="o",edgecolors="none")   
        if plot_only_coloc==False:
            im=ax.scatter(df.x1,df.y1,s=size,c=matrix[Gene1],cmap=c1,alpha=alpha,marker="o",edgecolors="none")
            im2=ax.scatter(df.x2,df.y2,s=size,c=matrix[Gene2],cmap=c2,alpha=alpha,marker="o",edgecolors="none")
        im3=ax.scatter(df.x4,df.y4,s=size,c=coloc,cmap=c3,alpha=alpha,marker="o",edgecolors="none")
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])


        if plot_only_coloc==False:
            axins = inset_axes(ax,
                       width=label_width,  # width = 5% of parent_bbox width
                       height=label_height,  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
            axins2 = inset_axes(ax,
                       width=label_width,  # width = 5% of parent_bbox width
                       height=label_height,  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0.15, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
            axins3 = inset_axes(ax,
                       width=label_width,  # width = 5% of parent_bbox width
                       height=label_height,  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0.3, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
            cb1=fig.colorbar(im, cax=axins,orientation='horizontal',shrink=0.3)
            cb2=fig.colorbar(im2, cax=axins2,orientation='horizontal',shrink=0.3)
            cb1.ax.set_title(label=Gene1,y=0.6,fontsize=fontsize)
            cb2.ax.set_title(label=Gene2,y=0.6,fontsize=fontsize)
            cb1.ax.tick_params(labelsize=fontsize)
            cb2.ax.tick_params(labelsize=fontsize)
           
        if plot_only_coloc==True:
            axins3 = inset_axes(ax,
                       width=label_width,  # width = 5% of parent_bbox width
                       height=label_height,  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.05, 0, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )

        textstr1=(str("%.2f" % per_total)+"% of Total")
        textstr2=(str("%.2f" % per_exp)+"% of Expressed")
        textstr3=(str("%.2f" % per_gene1)+"% of "+Gene1)
        textstr4=(str("%.2f" % per_gene2)+"% of "+Gene2)
        plt.text(1.0, 0.87, textstr1, fontsize=fontsize, transform=plt.gcf().transFigure)
        plt.text(1.0, 0.82, textstr2, fontsize=fontsize, transform=plt.gcf().transFigure)
        plt.text(1.0, 0.77, textstr3, fontsize=fontsize, transform=plt.gcf().transFigure)
        plt.text(1.0, 0.72, textstr4, fontsize=fontsize, transform=plt.gcf().transFigure)
        cb3=fig.colorbar(im3, cax=axins3,orientation='horizontal',shrink=0.3)
        cb3.ax.set_title(label="Mean Co-loc",y=0.6,fontsize=fontsize)
        cb3.ax.tick_params(labelsize=fontsize)
        del matrix["Gene1"]
        del matrix["Gene2"]
    def coloc_table(adata=adata,Gene1=Gene1,Gene2=Gene2,threshold=threshold):
            df=pd.DataFrame(adata.obsm["X_umap"],columns=["x","y"],index=adata.obs.index)
            df["x1"]=df["x"]
            df["x2"]=df["x"]
            df["x3"]=df["x"]
            df["x4"]=df["x"]
            matrix["Gene1"]=matrix[Gene1]
            matrix["Gene2"]=matrix[Gene2]
            coloc=(matrix["Gene1"]+matrix["Gene2"])/2
            df.x1[matrix.Gene1<threshold]="NA"
            df.x1[matrix.Gene2>=threshold]="NA"
            df.x2[matrix.Gene2<threshold]="NA"
            df.x2[matrix.Gene1>=threshold]="NA"
            df.x3[matrix.Gene1>=threshold]="NA"
            df.x3[matrix.Gene2>=threshold]="NA"         
            df.x4[matrix.Gene1<threshold]="NA"
            df.x4[matrix.Gene2<threshold]="NA"
            df.x1=df.x1.replace("NA",np.NaN)
            df.x2=df.x2.replace("NA",np.NaN)
            df.x3=df.x3.replace("NA",np.NaN)
            df.x4=df.x4.replace("NA",np.NaN)
            per_exp=df.x4.count()*100/(df.x1.count()+df.x2.count()+df.x4.count())
            per_total=df.x4.count()*100/(df.x1.count()+df.x2.count()+df.x4.count()+df.x3.count())
            per_gene1=df.x4.count()*100/(df.x1.count()+df.x4.count())
            per_gene2=df.x4.count()*100/(df.x2.count()+df.x4.count())
            name=Gene1+"_"+Gene2
            data={"Pair":name,
                  "Gene1":Gene1,
                  "Gene2":Gene2,
                  "Percentage_All":per_total,
                  "Percentage_Expressing":per_exp,
                  "Percentage_of_Gene1":per_gene1,
                  "Percentage_of_Gene2":per_gene2,
                  "Total_Cells":(df.x1.count()+df.x2.count()+df.x4.count()+df.x3.count()),
                  "Expressing_Cells":(df.x1.count()+df.x2.count()+df.x4.count()),
                  "Gene1_Cells":(df.x1.count()+df.x4.count()),
                  "Gene2_Cells":(df.x2.count()+df.x4.count())}
            del matrix["Gene1"]
            del matrix["Gene2"]
            return data

    if mode=="table":
        groups=[]
        pair=[]
        gene1=[]
        gene2=[]
        perall=[]
        perexp=[]
        perg1=[]
        perg2=[]
        tcells=[]
        expcells=[]
        g1cells=[]
        g2cells=[]
        if groupby!=None:
            group=set(list(adata.obs[groupby]))
            grouplist=[]
            for m in group:
                sub=adata[adata.obs[groupby]==m]
                for i,j in [(i,j) for i in Gene1 for j in Gene2]:
                    if (i!=j) and (f"{j}_{i}" not in groups):
                        x = coloc_table(sub,Gene1=i,Gene2=j)
                        grouplist=np.append(grouplist,m)
                        pair=np.append(pair,x["Pair"])
                        gene1=np.append(gene1,x["Gene1"])
                        gene2=np.append(gene2,x["Gene2"])
                        perall=np.append(perall,x["Percentage_All"])
                        perexp=np.append(perexp,x["Percentage_Expressing"])
                        perg1=np.append(perg1,x["Percentage_of_Gene1"])
                        perg2=np.append(perg2,x["Percentage_of_Gene2"])
                        tcells=np.append(tcells,x["Total_Cells"])
                        expcells=np.append(expcells,x["Expressing_Cells"])
                        g1cells=np.append(g1cells,x["Gene1_Cells"])
                        g2cells=np.append(g2cells,x["Gene2_Cells"])    
                        groups.append(f"{i}_{j}")
                        
            d=pd.DataFrame({
                    "Pair":pair,
                    "Gene1":gene1,
                    "Gene2":gene2,
                    "%_of_All":perall,
                    "%_of_Expressing":perexp,
                    "%_of_Gene1":perg1,
                    "%_of_Gene2":perg2,
                    "Total_Cells":tcells,
                    "Expressing_Cells":expcells,
                    "Gene1_Cells":g1cells,
                    "Gene2_Cells":g2cells,
                    groupby:grouplist})         
            d.reset_index(drop=True, inplace=True)
        else:
            for i,j in [(i,j) for i in Gene1 for j in Gene2]:
                if (i!=j) and (f"{j}_{i}" not in groups):
                    x = coloc_table(Gene1=i,Gene2=j)
                    pair=np.append(pair,x["Pair"])
                    gene1=np.append(gene1,x["Gene1"])
                    gene2=np.append(gene2,x["Gene2"])
                    perall=np.append(perall,x["Percentage_All"])
                    perexp=np.append(perexp,x["Percentage_Expressing"])
                    perg1=np.append(perg1,x["Percentage_of_Gene1"])
                    perg2=np.append(perg2,x["Percentage_of_Gene2"])
                    tcells=np.append(tcells,x["Total_Cells"])
                    expcells=np.append(expcells,x["Expressing_Cells"])
                    g1cells=np.append(g1cells,x["Gene1_Cells"])
                    g2cells=np.append(g2cells,x["Gene2_Cells"])    
                    groups.append(f"{i}_{j}")
            d=pd.DataFrame({
                    "Pair":pair,
                    "Gene1":gene1,
                    "Gene2":gene2,
                    "%_of_All":perall,
                    "%_of_Expressing":perexp,
                    "%_of_Gene1":perg1,
                    "%_of_Gene2":perg2,
                    "Total_Cells":tcells,
                    "Expressing_Cells":expcells,
                    "Gene1_Cells":g1cells,
                    "Gene2_Cells":g2cells
                   })  
            d.reset_index(drop=True, inplace=True)
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("End Time =", current_time)
        warnings.filterwarnings("ignore")   
        return d
    
    if mode=="plot":
        plot_coloc()
        
    if mode=="find_cells":
        df=pd.DataFrame(adata.obsm["X_umap"],columns=["x","y"],index=adata.obs.index)
        df["x"]="Yes"
        del df["y"]
        matrix["Gene1"]=matrix[Gene1]
        matrix["Gene2"]=matrix[Gene2]
        df.x[matrix.Gene1<threshold]="No"
        df.x[matrix.Gene2<threshold]="No"
        group=Gene1+"_"+Gene2+"_"+"Coloc"
        adata.obs[group]=df["x"]
        del df
        
    
sc.pp.filter_genes(adata, min_cells=25) # filter genes expressed in small number of cells

# filter genes  with low epxression
sc.pp.filter_genes(adata, min_counts=50) 
#previous filter means that at least 25 cells express the gene. 
#Setting min counts to 50 means each of those has to have at least 2 counts

adata.raw=adata

#create matrix after filtering
matrix=pd.DataFrame(adata.raw.X.toarray(),columns=adata.raw.var.index,index=adata.obs.index)
all1=np.asarray(adata.var.index)
all1=all1[:-1]
all2=np.asarray(adata.var.index)
df=gene_coloc(adata,mode="table",Gene1=all1,Gene2=all2,threshold=0.1)
df.to_csv(path_out+"dcm_ec_gene_coloc_all2.csv")