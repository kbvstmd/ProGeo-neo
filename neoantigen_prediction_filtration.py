import sys
import os
def handle_mkdir( ):
    A='mkdir  annovar_geneout'
    print A
    os.system(A)
def handle0(x): 
    O='perl'+' '+'soft/annovar/table_annovar.pl'+' '+x+' '+ 'soft/annovar/humandb/ -buildver hg38 -out annovar_geneout -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput'
    print O
    os.system(O)
def handle_mkdir1( ):
    B='mkdir  outfile2'
    print B
    os.system(B)    
def handle1( ):
    input1=open("annovar_geneout.hg38_multianno.txt",'r')
    output1=open("outfile2/gene-nonsynonymous-SNV.txt","w")
    for line in input1:
        ls1=line.strip().split("\t")
        if ls1[8]=='nonsynonymous SNV':
        
            ls2=ls1[9].strip('"').split(',')
            if len(ls2) >=2:
                for i in range(len(ls2)):
                    output1.write(ls2[i]+'\n')
                    #ls3=ls2[i].split(':')
            else:
                output1.write(ls1[9]+'\n')
def handle2( ):
    input1=open("reference_files/reuniprot.fasta","r")
    input2=open("outfile2/gene-nonsynonymous-SNV.txt","r")
    output1=open("outfile2/match1_proSeq.txt","w")
    dict1={}
    for line in input1:
        ls1=line.strip().split("\t")
        dict1[ls1[0]]=ls1[2]
        dict1[ls1[1]]=ls1[2]
    for line1 in input2:
        ls2=line1.strip().split(":")
        if dict1.has_key(ls2[0]) or dict1.has_key(ls2[1]):
            output1.write(ls2[0]+"\t"+ls2[4]+"\t"+dict1[ls2[0]]+"\n")
        else:
            pass
def handle_del( ):  
    lines_seen=set()
    outfile = open("outfile2/match_proSeq.txt","w")
    for line in open("outfile2/match1_proSeq.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
def handle_seq( ):
    input1=open("outfile2/match_proSeq.txt","r")
    output1=open("outfile2/gene-Var-proSeq.txt","w")
    for line in input1:
        ls1=line.strip().split("\t")
        ls2=ls1[1].split('.')
        ls3=ls2[1][0]
        ls4=ls2[1][-1]
        ls5=ls2[1][1:-1]
        ls6=int(ls5)-1
        ls7=list(ls1[2])
        if len(ls1[2]) >= int(ls2[1][1:-1]):  
            if ls7[ls6] == ls2[1][0]:
                ls7[ls6]=ls2[1][-1]
                seq="".join(ls7)
                output1.write(ls1[0]+'\t'+ls1[1]+'\t'+seq+'\n')
def handle3( ): 
    lines_seen=set()
    outfile = open("outfile2/Gene-Varseqence.csv","w")
    for line in open("outfile2/gene-Var-proSeq.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
def handle4( ): 
    input1=open("outfile2/Gene-Varseqence.csv","r")
    output1=open("outfile2/gene-21AA.fasta","w")
    output2=open("outfile2/21-AA.csv","w")
    for line in input1:
        ls=line.strip().split("\t")
        ls1=ls[1].split('.')
        pep=ls1[1][1:-1]
        seq=list(ls[2])
        pep1=seq[(int(pep)-11):(int(pep)+10)] 
        if len(pep1)==21:    
            seq1=''.join(pep1)
            output1.write('>'+ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+str(11)+"\n"+seq1+'\n') 
            output2.write(ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+str(11)+"\t"+seq1+'\n') 
        else:
            if int(pep)-11 <= 0:
                pep2=seq[0:(int(pep)+10)] 
                seq2=''.join(pep2)
                num1=pep
                output1.write('>'+ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+num1+"\n"+seq2+'\n')
                output2.write(ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+num1+"\t"+seq2+'\n')
            if  int(pep)+10 >= int(len(ls[2])):
                pep3=seq[(int(pep)-11):(int(pep)+1)]
                seq3=''.join(pep3)
                num2=str(11)
                output1.write('>'+ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+num2+"\n"+seq3+'\n')
                output2.write(ls[0]+':'+ls1[1][0]+ls1[1][-1]+':'+num2+"\t"+seq3+'\n')
def handle5(x): 
    w='netMHCpan -a '+x+' '+'outfile2/gene-21AA.fasta' + " > " + 'outfile2/my.out'
    print w
    os.system(w)
def handle6( ):
    input1=open("outfile2/my.out","r")
    output1=open("outfile2/my.csv","w")
    for line in input1:
        ls1=line.strip(" ").split('\00')
        ls2=ls1[0].strip(' ').split(' ')
        if "WB" in ls2[-1] or 'SB' in ls2[-1]:
            ls3='\t'.join(ls2[:])
            output1.write(ls3)
def handle7( ):
    import pandas as pd 
    infile ="outfile2/my.csv"
    outfile = "outfile2/myresult.csv"
    data = pd.read_table(infile, sep="\s+", header=None)
    temp = data.iloc[:, [1, 2,-5,-3,-1]]
    #print data.iloc[:, [-5]]
    temp.to_csv(outfile, sep="	", header=['HLA','Peptide','Gene','%Rank','BindLevel'], index=None)

def handle8( ):    
    input2=open("outfile2/21-AA.csv","r")  
    output=open("outfile2/myresult-noantiegn.csv",'w') 
    dict1={}
    for line in input2:
        ls1=line.strip('\n\t').split('\t')
        ls4=ls1[0].split(':')
        input1=open("outfile2/myresult.csv",'r')
        for line1 in input1:
            ls2=line1.strip().split('\t')
            if ls2[1] in ls1[1]:
                ls3=ls1[1].index(ls2[1])+1
                output.write(line1. strip()+'\t'+ls1[1]+'\t'+str(ls3)+'\t'+ls4[2]+'\n')               
            #print ls1[1]
              
def handle9( ):
    input1=open("outfile2/myresult-noantiegn.csv",'r')             
    output1=open("outfile2/candid1-noantiegn.txt",'w')             
    for line in input1:
        ls1=line.strip().split('\t')
        ls3=ls1[2].split('_')
        if float(ls1[6]) <= float(ls1[7]):
            if float(ls1[6])+float(len(ls1[1]))-1 >= float(ls1[7]):
                ls2=float(ls1[7])-float(ls1[6])+1
                #print ls2
                output1.write(ls1[0]+'\t'+ls3[0]+'\t'+ls1[1]+'\t'+ls1[3]+'\t'+ls1[4]+'\t'+ls3[1]+'\t'+str(ls2)+'\n')
def handle10( ): 
    lines_seen=set()
    outfile = open("candid-noantiegn.txt","w")
    for line in open("outfile2/candid1-noantiegn.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
def expression(x):
    index='soft/kallisto/kallisto  index  -i  outfile2/transcripts.idx  ' +' '+ x        
    print index
    os.system(index) 
def get_tpm(x,y):
    tpm='soft/kallisto/kallisto quant  -i  outfile2/transcripts.idx  -o  kallisto_output  -b  100 '  +x +' '+ y            
    print tpm
    os.system(tpm) 
def get_ID( ): 
    input1=open("reference_files/TRA_ENST_ID.txt","r")              
    input2=open("kallisto_output/abundance.tsv","r")                 
    output=open("outfile2/kallisto_gene_expression.txt","w")
    dic={}
    for line in input1:
        ls1=line.strip('\n').split('\t')
        dic[ls1[0]]=ls1[1]
    for line1 in input2:
        ls2=line1.strip('\n').split('\t')
        ls3=ls2[0].split('.')
        if ls2[4] > str(0):
            if dic.has_key(ls3[0]):
                output.write(dic[ls3[0]]+'\n')
def filter1( ): 
    input1=open("outfile2/kallisto_gene_expression.txt","r")   
    input2=open("candid-noantiegn.txt","r")                             
    output=open("candid_tpm_filter1-neoantiegn.txt","w")
    dic={}
    for line in input1:
        ls1=line.strip('\n').split('\t')
        dic[ls1[0]]= ' ' 
        for line1 in input2:
            ls2=line1.strip('\n').split('\t')
            if dic.has_key(ls2[1]):
                output.write(line1.strip('\n')+'\n')

# encoding=utf8
import io
import os
import click
import sys
#@click.command()
#@click.option('--raw_file_dir', '-raw', help='input raw data file')
#@click.option('--refer_file', '-ref', help='refer_file')
def main(raw_file_dir, refer_file):
    path=os.path.split(os.path.realpath(__file__))[0]  
    f = open(path+"/mqpar-lable.xml", "r")                       
    f1 = open(path+"/mqpar.xml", "w+")
    lines = f.readlines()
    reffasta = refer_file
    reffasta = "      <string>"+reffasta+"</string>\n"
    for i in range(0, 17):
        f1.write(lines[i])
    file_dir = raw_file_dir
    for filename in os.listdir(file_dir):
        if ".raw" in filename:
            tmp = "      <string>"+file_dir+"/"+filename+"</string>\n"
            print(tmp)
            f1.write(tmp)
    for i in range(18, 140):
        f1.write(lines[i])
    f1.write(reffasta)
    for i in range(141, 252):
        f1.write(lines[i])
    f1.close()
    cmd = "mono "+path+"/soft/MaxQuant/bin/MaxQuantCmd.exe "+path+"/mqpar.xml"
    print(cmd)
    os.system(cmd)
    return(raw_file_dir)
def handle11(x):
    input=open(x,'r')
    output=open('outfile2/Maxpep.txt','w')
    for line in input:
       ls1=line.strip('\n').split('\t')
       if ls1[0]=='Sequence':
           output.write('Sequence'+'\t'+'Gene'+'\t'+'AA Change'+'\t'+'Mut Position'+'\n')
       ls2=ls1[34].split(';')
       for i in range(len(ls2)):
           ls3=ls2[i].split('|')
           if len(ls3) >=2:
               ls4=ls3[1].split('.')
               ls5=ls4[1]
               ls6=ls5[-1]
               ls7=ls5[1:-1]
               #print ls7
               if float(ls1[36])<= int(str(ls7)) and int(str(ls7))<= float(ls1[37]):
                   ls8=float(ls7)-float(ls1[36])
                   if ls1[0][int(ls8)]==ls6:
                       ls9=ls5[0]+str(int(ls8)+1)+ls5[-1]
                        #output.write(ls1[0]+"\n")
                       output.write(ls1[0]+'\t'+ls3[0]+'\t'+ls9+'\t'+str(ls8+1)+'\t'+ls3[1]+'\n')
def handle12( ):
    input2=open('candid_tpm_filter1-noantiegn.txt',"r")
    output=open("outfile2/candid1_max_filter2-noantiegn.txt",'w') 
    for line in input2:
        ls1=line.strip('\n').split('\t')
        input1=open('outfile2/Maxpep.txt','r')
        for line1 in input1:
            ls2=line1.strip().split('\t')
            if ls1[2] in ls2[0]:
                output.write(line.strip()+'\n')
        input1.close()
    input2.close()
    output.close()
def handle13( ):
    lines_seen=set()
    outfile = open("candid_max_filter2-noantiegn.txt","w")
    for line in open("outfile2/candid1_max_filter2-noantiegn.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
def handle14( ):
    input1=open("candid_max_filter2-noantiegn.txt",'r')
    output1=open("candid_tpm_max_noantiegn.fasta",'w')
    for line in input1:
        ls1=line.strip('\n').split('\t')
        output1.write('>'+ls1[0]+'|'+ls1[1]+'|'+ls1[5]+'|'+ls1[6]+'\n'+ls1[2]+'\n')  
 
def blast( ):
    blast='soft/blast/bin/blastp  -query  candid_tpm_max_noantiegn.fasta   -db  iedb_db  -outfmt   "6 qacc qseq sacc sseq evalue length pident"  -evalue  100000000   -gapopen  11  -gapextend 1  >  candid_tpm_max_noantiegn-blast_out.txt'
    print blast
    os.system(blast) 
if __name__ == '__main__':
    temp=sys.argv                        
    handle_mkdir( )
    handle0(temp[1])  
    handle_mkdir1( )
    handle1( )     
    handle2( )
    handle_del( )
    handle_seq( )
    handle3( )
    handle4( )
    handle5(temp[2])   
    handle6( )
    handle7( )
    handle8( )
    handle9( )
    handle10( )
    expression(temp[3])
    get_tpm(temp[4],temp[5])
    get_ID( )
    filter1( )
    re_path=main(temp[6],temp[7])
    print(re_path)	
    handle11(re_path+'/combined/txt/peptides.txt')
    handle12( )
    handle13( )
    handle14( )
    blast( )
