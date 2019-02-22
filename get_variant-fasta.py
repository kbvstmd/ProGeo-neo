import sys
import os
def handle0( ):
    Z='mkdir outfile1'
    print Z
    os.system(Z)   
def handle1( ):
    A='soft/bwa/bwa index reference_files/hg38.fasta'
    print A
    os.system(A)
def handle2(y,z):
    B='soft/bwa/bwa mem -t 8 -R "@RG\\tID:foo\\tSM:bar\\tLB:library1" reference_files/hg38.fasta'+' '+y+' '+z+' '+'>'+' outfile1/sample.sam'
    print B
    os.system(B)
    #handle2('fastq1','fastq2 )
def handle3( ):
    C='soft/samtools/samtools fixmate -O bam outfile1/sample.sam outfile1/sample.bam'  
    print C
    os.system(C)
    #handle3('', '')
def handle4( ):
    D='mkdir tmp'
    print D
    os.system(D)
def handle5( ):
    E='soft/samtools/samtools sort -O bam -o outfile1/sorted_sample.bam -T  temp outfile1/sample.bam' 
    os.system(E)
    print E
    #handle5('sorted_sample.bam','sample.bam')
def handle6( ):   
    F='soft/gatk/gatk MarkDuplicates -I outfile1/sorted_sample.bam -O outfile1/sorted_markedup.bam -M outfile1/sorted_markedup_metric.bam'
    print F
    os.system(F)
    #handle6('sorted-sample.bam', 'sorted_markedup.bam', 'sorted_markedup_metric.bam')
def handle7( ):
    G='soft/gatk/gatk AddOrReplaceReadGroups -I  outfile1/sorted_markedup.bam -O  outfile1/sorted_markedup_Add.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM Cancer' 
    print G
    os.system(G)    
    #handle7('sorted_markedup.bam', 'sorted_markedup_Add.bam')
def handle8( ):
    H='soft/gatk/gatk BaseRecalibrator -I outfile1/sorted_markedup_Add.bam  -R  reference_files/hg38.fasta  --known-sites  reference_files/dbsnp_146.hg38.vcf -O outfile1/recal_data-sample.table'  
    print H
    os.system(H)
    #handle8('sorted_markedup_Add.bam','reference_files/hg38.fasta',"dbsnp_146.hg38.vcf",'recal_data-sample.table')
def handle9( ):
    I='soft/gatk/gatk ApplyBQSR -R  reference_files/hg38.fasta -I  outfile1/sorted_markedup_Add.bam --bqsr-recal-file  outfile1/recal_data-sample.table -O  outfile1/recal-sample.bam'
    print I
    os.system(I)
    #handle9('reference_files/hg38.fasta','sorted_markedup_Add.bam','recal_data-sample.table','recal-sample.bam')
def handle10( ): 
    J='soft/samtools/samtools index  outfile1/recal-sample.bam '
    print J
    os.system(J)
def handle11( ):
    K='soft/bcftools/bcftools mpileup -Ou -f reference_files/hg38.fasta outfile1/recal-sample.bam | soft/bcftools/bcftools call -vmO z -o  outfile1/sample.vcf.gz'
    print K
    os.system(K)
    #handle11('reference_files/hg38.fasta','recal-sample.bam','sample.vcf.gz')
def handle12( ):
    L='soft/bcftools/bcftools tabix -p vcf outfile1/sample.vcf.gz'
    print L
    os.system(L)
    #handle12('sample.vcf.gz')
def handle13( ):
    M="soft/bcftools/bcftools filter -O z -o outfile1/sample.filtered.vcf.gz -s LOWQUAL -i '%QUAL>20' outfile1/sample.vcf.gz"
    print M
    os.system(M)
    #handle13('sample.filtered.vcf.gz','sample.vcf.gz')
 
def handle14( ):
    Q='gunzip outfile1/sample.filtered.vcf.gz'
    os.system(Q)
#handle14(x)
def handle15( ):
    N='mkdir annovar_out'
    print N
    os.system(N)
def handle16( ): 
    O='perl'+' '+'soft/annovar/table_annovar.pl'+' '+'outfile1/sample.filtered.vcf soft/annovar/humandb/ -buildver hg38 -out annovar_out -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput'
    print O
    os.system(O)
    #handle15('sample2.filtered.vcf.gz','annovar_tmp')
def handle17( ):
    input1=open('annovar_out.hg38_multianno.txt','r')
    output1=open("outfile1/sample-nonsynonymous-SNV.txt","w")
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
def handle18( ):
    input1=open('reference_files/uniprot.fasta','r')
    output=open('outfile1/reuniprot.fasta','w')
    a=''
    b=''
    for line in input1:    
        if '>' in line:
            if a !='':
                output.write(a+'\n')      
            b=line.split('|')[2].split('_')[0]  
            c= line.split(' ')
            d=c[-3].split('=')
        #if 
            e=d[1]
        #print e 
            output.write(b+'\t'+e+'\t')
            a=''
        else:
            a+=line.strip()
    output.write(a+'\n')
def handle19( ):
    input1=open("outfile1/reuniprot.fasta","r")
    input2=open("outfile1/sample-nonsynonymous-SNV.txt","r")
    output1=open("outfile1/candidate_match_proSeq.txt","w")
    dict1={}
    dict2={}
    for line in input1:
        ls1=line.strip().split("\t")
        dict1[ls1[0]]=ls1[2]
        dict2[ls1[1]]=ls1[2]
    #print ls1[0].strip(">")
    #print ls1[1]
    for line1 in input2:
        ls2=line1.strip().split(":")
    
        if dict1.has_key(ls2[0]) or dict1.has_key(ls2[1]):
            output1.write(ls2[0]+"\t"+ls2[4]+"\t"+dict1[ls2[0]]+"\n")
        else:
            pass
def handle20( ):  
    lines_seen=set()
    outfile = open("outfile1/match_proSeq.txt","w")
    for line in open("outfile1/candidate_match_proSeq.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
    
def handle21( ):
    input1=open("outfile1/match_proSeq.txt","r")
    output1=open("Var-proSeq.fasta","w")
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
                output1.write(">"+ls1[0]+'|'+ls1[1]+'\n'+seq+'\n')

if __name__ == '__main__':
    temp=sys.argv   
    handle0( )
    handle1( ) 
    handle2(temp[1],temp[2])
    handle3( )
    handle4( )
    handle5( )
    handle6( )
    handle7( )
    handle8( )
    handle9( )
    handle10( )
    handle11( )
    handle12( )
    handle13( )
    handle14( )
    handle15( )
    handle16( )
    handle17( )
    handle18( )
    handle19( )
    handle20( )
    handle21( )
