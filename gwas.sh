#注：该流程目前不适用于非人物种
#注：该流程目前不适合于家系类样本的GWAS分析
#注：该流程只能用于散发无亲缘关系类样本的GWAS数据分析
#注：该流程的impute数据必须是Genomestudio导出的top allele的结果。或者是Genomestudio安装的plink插件导出的结果。否则imputation可能会出错。
#注：该流程运行前必须将所有的case和control样本整理好,生成二进制格式的输入文件.
#注：整理输入文件时可能用到如下命令，请学会plink软件的最基本的操作。
#注：--make-bed --recode --vcf --allow-no-sex --merge --bmerge --update-sex --update-map --update-name --flip --keep --remove --extrect --exclude --filter-cases --filter-controls等，请熟悉这些命令的详细用途。
#注：FID和IID必须一样，且唯一。
workdir=/lustre/public_data/SNP/beiyi/final/
shapeit=/opt/shapeit.v2.r837.GLIBCv2.12.Linux.static/bin/shapeit
gcta=/opt/gcta_1.26.0/gcta64
annovar=/opt/annovar
kg=/lustre/reference/1KG_phase3
inpute=beiyi     #需要提供plink的inpute文件，如beiyi.bed beiyi.bim beiyi.fam.
maf=0.01
geno=0.05
mind=0.05
missingP=0.001
comm1=--allow-no-sex
hwe=0.001
threads=10
sex_check=0     #1 or 0 #是否需要将sex有问题的样本进行去除。1为去除，即默认，0为不去除。
Pi_HAT=0.15      #0.5表明样本是同胞对时的临界值，这里取到0.4
LD_region=/lustre/home/zhangyongbiao/bin/highLDregions4bim_b37.awk
strand_file=/lustre/public_data/SNP/beiyi/strand_file/GSA-24v1-0_A1-b37.strand
impute_command=/lustre/reference/1KG_phase3/human_b37.txt
snptest=/opt/snptest_v2.5.4-beta3_linux_x86_64_dynamic/snptest_v2.5.4-beta3
minimum_predictor_count=120      #该值=sample_number*2*0.01
subgroup="ALL A B H"          #提供所有phenotype的name，注意必须和pheno or sample文件中的名字一致,注意，ALL即所有样本也是其一子集。
rmind=rmind.txt				#分析开始之前就要删除的那些无表型、项目不相关等样本
pheno_file=/lustre/public_data/SNP/beiyi/final/pheno.txt #这里列出客户提供的所有协变量，其中必须包含age和sex（注意title一律用小写）。


#如下参数是针对北医项目特有的，因为北医的imputation时，case和control是分开impute的
#case_gen_dir=/lustre/public_data/SNP/beiyi/imputation/case/imputation
control_gen_dir=/lustre/public_data/SNP/beiyi/imputation/control/imputation
case_pheno=/lustre/public_data/SNP/beiyi/final/case.sample		#该文件需要手动生成
control_pheno=/lustre/public_data/SNP/beiyi/final/control.sample		#该文件需要手动生成
exclude_sample=/lustre/public_data/SNP/beiyi/final/exclude_sample

#当存在多个分组时，用整体计算出来的pca并不能满足亚族人群分层校正时获得lambda值，因此需要对每个亚组+control的数据重新进行pca的分析。
A=/lustre/public_data/SNP/beiyi/final/A
B=/lustre/public_data/SNP/beiyi/final/B
H=/lustre/public_data/SNP/beiyi/final/H
A_control=/lustre/public_data/SNP/beiyi/final/A_control
B_control=/lustre/public_data/SNP/beiyi/final/B_control
H_control=/lustre/public_data/SNP/beiyi/final/H_control
ALL_control=/lustre/public_data/SNP/beiyi/final/ALL_control
#每个pheno的prevalence rate
A_prevalence=
B_prevalence=
H_prevalence=
ALL_prevalence=



##数据分析前质控

##基本质控
mkdir genotype
mv $inpute* genotype
cd genotype

echo "
###去除项目那些不相关、无表型等样本
plink --bfile $inpute --remove ../$rmind --make-bed --out $inpute.rmind $comm1
##去除maf设定值的SNP位点
plink --bfile $inpute.rmind --maf $maf --make-bed --out $inpute.rmind.maf $comm1
##去除低于设定成功率的位点
plink --bfile $inpute.rmind.maf --geno $geno --make-bed --out $inpute.rmind.fam.geno $comm1
##去除在各亚组中成功率低于设定值的SNP位点
for subgroup in $subgroup; do
plink --bfile $inpute.rmind.fam.geno --geno $geno --keep ../$subgroup --make-bed --out $inpute.rmind.fam.geno.subgroup_${subgroup}
done

#cat *subgroup*bim | awk '{print $2}' | uniq - > geno.filter
#grep -F -f *subgroup*bim |awk '{print $2}' |uniq - > geno.filter

plink --bfile $inpute.rmind.fam.geno --extract $inpute.rmind.fam.geno.subgroup_A.bim --make-bed --out $inpute.rmind.fam.geno2 $comm1
plink --bfile $inpute.rmind.fam.geno2 --extract $inpute.rmind.fam.geno.subgroup_B.bim --make-bed --out $inpute.rmind.fam.geno3 $comm1
plink --bfile $inpute.rmind.fam.geno3 --extract $inpute.rmind.fam.geno.subgroup_H.bim --make-bed --out $inpute.rmind.fam.geno4 $comm1

"

echo "
##去除低于设定成功率的样本,生成clean数据
plink --bfile $inpute.rmind.fam.geno4 --mind $mind --make-bed --out clean $comm1
#单独生成case和control的数据集
plink --bfile clean --filter-cases --make-bed --out $inpute.cases
plink --bfile clean --filter-controls --make-bed --out $inpute.controls
"
cd ..



##其他质控
mkdir QC
cd QC

echo "
#性别质控
mkdir check_sex
plink --bfile ../genotype/clean --check-sex --out check_sex/clean    #sex_check
grep "PROBLEM" check_sex/clean.sexcheck | awk '{print $1,$2}' > check_sex/sex.error #提取有问题的样本
if [$sex_check -eq 1]
then
	plink --bfile ../genotype/clean --remove check_sex/sex.error --make-bed --out ../genotype/clean.sex  #去除性别有问题的样本
else
	plink --bfile ../genotype/clean --make-bed --out ../genotype/clean.sex
fi
"

echo "
#HWE质控
mkdir HWE
plink --bfile ../genotype/clean.sex --hardy --out HWE/$inpute $comm1
#剔除control中不符合HWE设定值的SNP位点
grep "UNAFF" HWE/beiyi.hwe | awk '$9 <= 0.001 {print $2}' > HWE/SNPs.departure.hwe
plink --bfile ../genotype/clean.sex --exclude HWE/SNPs.departure.hwe --make-bed --out ../genotype/clean.sex.hwe


#missing rate在case和control间差异显著性质控
mkdir missing
plink --bfile ../genotype/clean.sex.hwe --test-missing --out missing/clean.sex.hwe
awk '$5 <= 0.001 {print $2}' missing/clean.sex.hwe.missing > missing/missing.diff.snp
plink --bfile ../genotype/clean.sex.hwe --exclude missing/missing.diff.snp --make-bed --out ../genotype/clean.sex.hwe.missing $comm1
plink --bfile ../genotype/clean.sex.hwe.missing --make-bed --out ../genotype/clean.sex.hwe $comm1


#missing rate 在亚组合control间差异显著性控制
plink --bfile ../genotype/clean.sex.hwe --keep $A_control --test-missing --out missing/clean.sex.hwe
awk '$5 <= 0.001 {print $2}' missing/clean.sex.hwe.missing > missing/missing.diff.snp
plink --bfile ../genotype/clean.sex.hwe --exclude missing/missing.diff.snp --make-bed --out ../genotype/clean.sex.hwe.missing $comm1
plink --bfile ../genotype/clean.sex.hwe.missing --make-bed --out ../genotype/clean.sex.hwe $comm1

plink --bfile ../genotype/clean.sex.hwe --keep $B_control --test-missing --out missing/clean.sex.hwe
awk '$5 <= 0.001 {print $2}' missing/clean.sex.hwe.missing > missing/missing.diff.snp
plink --bfile ../genotype/clean.sex.hwe --exclude missing/missing.diff.snp --make-bed --out ../genotype/clean.sex.hwe.missing $comm1
plink --bfile ../genotype/clean.sex.hwe.missing --make-bed --out ../genotype/clean.sex.hwe $comm1

plink --bfile ../genotype/clean.sex.hwe --keep $H_control --test-missing --out missing/clean.sex.hwe
awk '$5 <= 0.001 {print $2}' missing/clean.sex.hwe.missing > missing/missing.diff.snp
plink --bfile ../genotype/clean.sex.hwe --exclude missing/missing.diff.snp --make-bed --out ../genotype/clean.sex.hwe.missing $comm1
plink --bfile ../genotype/clean.sex.hwe.missing --make-bed --out ../genotype/clean.sex.hwe $comm1
"

echo "
#亲缘关系质控
mkdir relation
#对于位于同一个LD内的SNP位点进行purn
plink --bfile ../genotype/clean.sex.hwe --indep-pairwise 50 5 0.2 --out ../genotype/clean.sex.hwe $comm1
plink --bfile ../genotype/clean.sex.hwe --extract ../genotype/clean.sex.hwe.prune.in --make-bed --out ../genotype/clean.sex.hwe.purned $comm1
#提取染色体中高度连锁的区域,并去除高度连锁不平衡区域内的位点
awk -f $LD_region ../genotype/clean.sex.hwe.purned.bim > relation/highLDexcludes
awk '($1 < 1) || ($1 > 22) {print $2}' ../genotype/clean.sex.hwe.purned.bim  > relation/autosomeexcludes
cat relation/highLDexcludes relation/autosomeexcludes > relation/highLD_and_autosomal_excludes
plink --bfile ../genotype/clean.sex.hwe.purned --exclude relation/highLD_and_autosomal_excludes --make-bed --out ../genotype/clean.sex.hwe.purned.rmhighLD $comm1

#计算Pi-hat并去除亲缘关系较近的样本
#plink --bfile ../genotype/clean.sex.hwe.purned.rmhighLD --genome --out ../genotype/clean.sex.hwe.purned.rmhighLD.IBD $comm1
awk '$10 >= 0.15 {print $0}' ../genotype/clean.sex.hwe.purned.rmhighLD.IBD.genome > relation/related_individuals
awk '{print $3," ",$4}' relation/related_individuals > relation/individual.violate.IBD
plink --bfile ../genotype/clean.sex.hwe.purned.rmhighLD --remove relation/individual.violate.IBD --make-bed --out ../genotype/clean.sex.hwe.purned.rmhighLD.rmrelated $comm1
plink --bfile ../genotype/clean.sex.hwe --remove relation/individual.violate.IBD --make-bed --out ../genotype/clean.sex.hwe.rmrelated
#"

echo "
###去除outlier
mkdir outlier
#####remove outliers with neighbour in plink. we first get indivuduals with Z score < -2, then pick indivudials have long distance to all 5 exomined individuals
plink --bfile ../genotype/clean.sex.hwe.purned.rmhighLD.rmrelated --cluster --neighbour 1 5 --out outlier/clean.sex.hwe.purned.rmhighLD.rmrelated
awk '$5 <= -2 {print $1," ",$2}' outlier/clean.sex.hwe.purned.rmhighLD.rmrelated.nearest | sort | uniq -cd | sort -nr | awk '$1 ==5 {print $2," ",$3}' > outlier/outliers
plink --bfile ../genotype/clean.sex.hwe.rmrelated --remove outlier/outliers --make-bed --out ../genotype/clean.sex.hwe.rmrelated.rmoutlier
"

##去除首次GWAS分析假阳性位点
#生成mds作为协变量,使用mds作为协变量进行关联分析,提取所有P<=1E10-7的SNP位点
mkdir mds
mkdir assoc_qc
echo "
for pheno in $subgroup; do
#plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --keep ../${pheno}_control --mds-plot 10 --cluster --out mds/clean.sex.hwe.rmrelated.rmoutlier.$pheno --threads 20
#plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --logistic --out assoc_qc/clean.sex.hwe.rmrelated.rmoutlier.$pheno --covar mds/clean.sex.hwe.rmrelated.rmoutlier.${pheno}.mds --covar-name C1-C10 --hide-covar --pheno $pheno_file --pheno-name $pheno --threads 20
awk '$9 <= 0.0000001 {print $2}' assoc_qc/clean.sex.hwe.rmrelated.rmoutlier.${pheno}.assoc.logistic > assoc_qc/${pheno}.sig
done

cat assoc_qc/*.sig |sort - | uniq - > ../manually.sheck.snps.txt
#手动检查所有这些显著的变异位点，并将那些假阳性位点手动剔出。
#
#plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --exclude ../false.postive.snps --make-bed --out ../genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP

#"
#生成针对每个亚族的pheno文件，即包含了age，sex，cov等
echo "
workdir = '$workdir'
subgroups = '${subgroup}'
subgroup=list(subgroups.strip().split())
print subgroup
mypheno= open (workdir + 'pheno.txt')
mydict={}
for line in mypheno:
        tokens=line.strip().split()
        if len(tokens) <=2:
                continue
        tempstr=''
        for i in range(2,len(tokens)):
                tempstr+=tokens[i]+'\t'
        mydict[tokens[0]]=tempstr
print len(mydict.keys())
for item in subgroup:
        myfile=open(workdir + 'QC/mds/clean.sex.hwe.rmrelated.rmoutlier.'+item+'.mds')
        myout = open (workdir + item+'.temp.pheno','w')
        for line in myfile:
                tokens=line.strip().split()
                if tokens[0] not in mydict.keys():
                        continue
		tempstr=''
		del tokens[2]
		for i in tokens:
			tempstr+=i+'\t'
                myout.write(tempstr+mydict[tokens[0]].strip()+'\n')
        myfile.close()
        myout.close()

" > ../generate.pheno.step1.py

echo "
workdir = '$workdir'
subgroups = '${subgroup}'
subgroup=list(subgroups.strip().split())
print subgroup
for item in subgroup:
	myfile = open(workdir+item+'.temp.pheno')
	myout = open(workdir+item+'.pheno','w')
	myout2 = open (workdir+item+'.plink.pheno','w')
	mydict={}	
	for line in myfile:
		tokens=line.strip().split()
		mydict[tokens[0]]=line.strip()
	myALL=open(workdir+'ALL.temp.pheno')
	for line in myALL:
		if 'FID' in line:
			myout2.write(line)
			line=line.replace('FID','ID_1')
			line=line.replace('IID','ID_2')
			myout.write(line+'0\t0\tC\tC\tC\tC\tC\tC\tC\tC\tC\tC\tC\tC\t'+len(subgroup)*'P\t'+'\n')
			continue
		tokens=line.strip().split()
		if tokens[0] in mydict.keys():
			myout.write(mydict[tokens[0]]+'\n')
			myout2.write(mydict[tokens[0]]+'\n')
		else:
			myout.write(line)
			myout2.write(line)
	myfile.close()
	myout.close()
"  > ../generate.pheno.step2.py

cd ..

#python generate.pheno.step1.py
#python generate.pheno.step2.py
rm *.temp.pheno

##############################################注意,这里不能一次运行结束,因为需要手动检查位点的可信度.同时要对生成的新的pheno文件按照SNPtest的格式进行调整，主要调整的就是第二行#######################################
##############################################第二行为第一行的解释行，通过单个字母解释且一一对应，如果是协变量，对应填写C，如果是bin数据，对应填写B，如果是pheno对应填写P，其他均为0#################################


#####imputation######
echo "
##imputation前准备
#负链到正链
#首先提取所有位于负链的SNP位点,然后flip，
awk '{if($5 == "-") print $1}' $strand_file > genotype/SNP.minus.strand
plink --bfile genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP --flip genotype/SNP.minus.strand --make-bed --out genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand
#去除相同物理位置的SNP位点
#"

echo "
workdir = '$workdir'
myfile = open (workdir + 'genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand.bim')
myout = open (workdir + 'genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand.new.bim','w')
mydict={}
for line in myfile:
        tokens = line.strip().split()
        if tokens[0] == '0':
                continue
        mydict[tokens[0]+':'+tokens[3]] = line
for item in mydict.keys():
        myout.write(mydict[item].strip().split()[1]+'\n')
myfile.close()
myout.close()
" > rm.same.pos.py

echo "
python rm.same.pos.py
plink --bfile genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand --extract genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand.new.bim --make-bed --out genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand.rmpos $comm1
#"

mkdir imputation
echo "
#分割染色体
for i  in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
plink --bfile genotype/clean.sex.hwe.rmrelated.rmoutlier.rmFP.plusstrand.rmpos --chr $i --make-bed --out imputation/clean.chr$i
done
#"
cd imputation
#使用shapeit开始进行单体型构建
mkdir shapeit
mkdir chunks
mkdir impute2

echo "
$shapeit -B clean.chr1 -M ${kg}/genetic_map_chr1_combined_b37.txt -O shapeit/clean.chr1.phased -T 30 &
$shapeit -B clean.chr2 -M ${kg}/genetic_map_chr2_combined_b37.txt -O shapeit/clean.chr2.phased -T 30 &
$shapeit -B clean.chr3 -M ${kg}/genetic_map_chr3_combined_b37.txt -O shapeit/clean.chr3.phased -T 30 &
$shapeit -B clean.chr4 -M ${kg}/genetic_map_chr4_combined_b37.txt -O shapeit/clean.chr4.phased -T 30 &
$shapeit -B clean.chr5 -M ${kg}/genetic_map_chr5_combined_b37.txt -O shapeit/clean.chr5.phased -T 30 &
$shapeit -B clean.chr6 -M ${kg}/genetic_map_chr6_combined_b37.txt -O shapeit/clean.chr6.phased -T 30 &
$shapeit -B clean.chr7 -M ${kg}/genetic_map_chr7_combined_b37.txt -O shapeit/clean.chr7.phased -T 30 &
$shapeit -B clean.chr8 -M ${kg}/genetic_map_chr8_combined_b37.txt -O shapeit/clean.chr8.phased -T 30 &
$shapeit -B clean.chr9 -M ${kg}/genetic_map_chr9_combined_b37.txt -O shapeit/clean.chr9.phased -T 30 &
$shapeit -B clean.chr10 -M ${kg}/genetic_map_chr10_combined_b37.txt -O shapeit/clean.chr10.phased -T 30 &
$shapeit -B clean.chr11 -M ${kg}/genetic_map_chr11_combined_b37.txt -O shapeit/clean.chr11.phased -T 30 &
$shapeit -B clean.chr12 -M ${kg}/genetic_map_chr12_combined_b37.txt -O shapeit/clean.chr12.phased -T 30 &
$shapeit -B clean.chr13 -M ${kg}/genetic_map_chr13_combined_b37.txt -O shapeit/clean.chr13.phased -T 30 &
$shapeit -B clean.chr14 -M ${kg}/genetic_map_chr14_combined_b37.txt -O shapeit/clean.chr14.phased -T 30 &
$shapeit -B clean.chr15 -M ${kg}/genetic_map_chr15_combined_b37.txt -O shapeit/clean.chr15.phased -T 30 &
$shapeit -B clean.chr16 -M ${kg}/genetic_map_chr16_combined_b37.txt -O shapeit/clean.chr16.phased -T 30 &
$shapeit -B clean.chr17 -M ${kg}/genetic_map_chr17_combined_b37.txt -O shapeit/clean.chr17.phased -T 30 &
$shapeit -B clean.chr18 -M ${kg}/genetic_map_chr18_combined_b37.txt -O shapeit/clean.chr18.phased -T 30 &
$shapeit -B clean.chr19 -M ${kg}/genetic_map_chr19_combined_b37.txt -O shapeit/clean.chr19.phased -T 30 &
$shapeit -B clean.chr20 -M ${kg}/genetic_map_chr20_combined_b37.txt -O shapeit/clean.chr20.phased -T 30 &
$shapeit -B clean.chr21 -M ${kg}/genetic_map_chr21_combined_b37.txt -O shapeit/clean.chr21.phased -T 30 &
$shapeit -B clean.chr22 -M ${kg}/genetic_map_chr22_combined_b37.txt -O shapeit/clean.chr22.phased -T 30 
"
echo "

#开始imputation,这里默认投放20个任务。如果要投放更多任务加速进程，比如40个，请把下方注释的#取出。
split -n l/40 -d $impute_command #如果要投放40个任务，请将-n l/20 改为 -n l/40 

sh x00 & sh x01 & sh x02 & sh x03 & sh x04 & sh x05 & sh x06 & sh x07 & sh x08 & sh x09 & sh x10 & sh x11 & sh x12 & sh x13 & sh x14 & sh x15 & sh x16 & sh x17 & sh x18 & sh x19 & sh x20 & sh x21 & sh x22 & sh x23 & sh x24 & sh x25 & sh x26 & sh x27 & sh x28 & sh x29 & sh x30 & sh x31 & sh x32 & sh x33 & sh x34 & sh x35 & sh x36 & sh x37 & sh x38 & sh x39 #如果要投放40个任务，请将该行的第一个#删除。
#$imputation_run
"
echo "
#合并gz文件
zcat chunks/clean.chr1.chunk*.gz | gzip > impute2/clean.chr1.gz &
zcat chunks/clean.chr2.chunk*.gz | gzip > impute2/clean.chr2.gz &
zcat chunks/clean.chr3.chunk*.gz | gzip > impute2/clean.chr3.gz &
zcat chunks/clean.chr4.chunk*.gz | gzip > impute2/clean.chr4.gz &
zcat chunks/clean.chr5.chunk*.gz | gzip > impute2/clean.chr5.gz &
zcat chunks/clean.chr6.chunk*.gz | gzip > impute2/clean.chr6.gz &
zcat chunks/clean.chr7.chunk*.gz | gzip > impute2/clean.chr7.gz &
zcat chunks/clean.chr8.chunk*.gz | gzip > impute2/clean.chr8.gz &
zcat chunks/clean.chr9.chunk*.gz | gzip > impute2/clean.chr9.gz &
zcat chunks/clean.chr10.chunk*.gz | gzip > impute2/clean.chr10.gz &
zcat chunks/clean.chr11.chunk*.gz | gzip > impute2/clean.chr11.gz &
zcat chunks/clean.chr12.chunk*.gz | gzip > impute2/clean.chr12.gz &
zcat chunks/clean.chr13.chunk*.gz | gzip > impute2/clean.chr13.gz &
zcat chunks/clean.chr14.chunk*.gz | gzip > impute2/clean.chr14.gz &
zcat chunks/clean.chr15.chunk*.gz | gzip > impute2/clean.chr15.gz &
zcat chunks/clean.chr16.chunk*.gz | gzip > impute2/clean.chr16.gz &
zcat chunks/clean.chr17.chunk*.gz | gzip > impute2/clean.chr17.gz &
zcat chunks/clean.chr18.chunk*.gz | gzip > impute2/clean.chr18.gz &
zcat chunks/clean.chr19.chunk*.gz | gzip > impute2/clean.chr19.gz &
zcat chunks/clean.chr20.chunk*.gz | gzip > impute2/clean.chr20.gz &
zcat chunks/clean.chr21.chunk*.gz | gzip > impute2/clean.chr21.gz &
zcat chunks/clean.chr22.chunk*.gz | gzip > impute2/clean.chr22.gz 
#"
cd ..

echo "

#SNPtest
echo "
cd imputation
#mkdir snptest

echo "
for name in $subgroup; do
$snptest -data impute2/clean.chr1.gz ../${name}.pheno -o snptest/clean.chr1.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr1.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr2.gz ../${name}.pheno -o snptest/clean.chr2.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr2.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr3.gz ../${name}.pheno -o snptest/clean.chr3.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr3.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr4.gz ../${name}.pheno -o snptest/clean.chr4.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr4.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr5.gz ../${name}.pheno -o snptest/clean.chr5.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr5.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr6.gz ../${name}.pheno -o snptest/clean.chr6.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr6.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr7.gz ../${name}.pheno -o snptest/clean.chr7.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr7.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr8.gz ../${name}.pheno -o snptest/clean.chr8.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr8.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr9.gz ../${name}.pheno -o snptest/clean.chr9.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr9.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr10.gz ../${name}.pheno -o snptest/clean.chr10.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr10.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr11.gz ../${name}.pheno -o snptest/clean.chr11.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr11.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr12.gz ../${name}.pheno -o snptest/clean.chr12.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr12.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr13.gz ../${name}.pheno -o snptest/clean.chr13.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr13.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr14.gz ../${name}.pheno -o snptest/clean.chr14.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr14.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr15.gz ../${name}.pheno -o snptest/clean.chr15.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr15.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr16.gz ../${name}.pheno -o snptest/clean.chr16.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr16.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr17.gz ../${name}.pheno -o snptest/clean.chr17.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr17.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr18.gz ../${name}.pheno -o snptest/clean.chr18.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr18.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr19.gz ../${name}.pheno -o snptest/clean.chr19.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr19.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr20.gz ../${name}.pheno -o snptest/clean.chr20.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr20.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr21.gz ../${name}.pheno -o snptest/clean.chr21.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr21.${name}_snptest.log -exclude_samples $exclude_sample  &
$snptest -data impute2/clean.chr22.gz ../${name}.pheno -o snptest/clean.chr22.${name}.snptest -frequentist 1 2 3 4 5 -method score -pheno ${name} -hwe -cov_names age sex C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 -log snptest/clean.chr22.${name}_snptest.log -exclude_samples $exclude_sample  &

done


#"







cd ..


#使用plink进行常规关联分析

echo "

mkdir association
cd association
mkdir allelic
mkdir model
mkdir logistic
mkdir LD

for pheno in ${subgroup}; do
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --assoc --out allelic/clean.sex.hwe.rmrelated.rmoutlier.$pheno --ci 0.95 --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --model trend --out model/clean.sex.hwe.rmrelated.rmoutlier.trend.$pheno --ci 0.95 --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --model trend-only --out model/clean.sex.hwe.rmrelated.rmoutlier.trend-only.$pheno --ci 0.95 --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --model dom --out model/clean.sex.hwe.rmrelated.rmoutlier.dom.$pheno --ci 0.95 --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --model rec --out model/clean.sex.hwe.rmrelated.rmoutlier.rec.$pheno --ci 0.95 --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --logistic --out logistic/clean.sex.hwe.rmrelated.rmoutlier.$pheno --ci 0.95 --covar ../${pheno}.plink.pheno --covar-name age,sex,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --hide-covar --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --exclude ../false.postive.snps
plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --logistic --out logistic/clean.sex.hwe.rmrelated.rmoutlier.beta.$pheno --ci 0.95 --beta --adjust --pheno ../${pheno}.plink.pheno --pheno-name $pheno --covar ../${pheno}.plink.pheno --covar-name age,sex,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --hide-covar --exclude ../false.postive.snps
done


plink --bfile ../genotype/clean.sex.hwe.rmrelated.rmoutlier --r2 --out LD/clean.sex.hwe.rmrelated.rmoutlier --ld-window-r2 0.2 --exclude ../false.postive.snps


#"

#使用tcga进行表型变异解释度分析
ehco "
mkdir gcta
for name in $subgroup;do
awk '{print $1,$2,$15}' ${name}.plink.pheno > gcta/clean_phenoforgcta_${name}.txt
done
#for name in $subgroup; do
	#make a grm for each chromosome using the unimputed bed file. The entire analysis should also be repeated using the imputed data
	for chr in $(seq 1 22); do 
	$gcta --bfile imputation/clean.chr$chr --make-grm --out gcta/clean_ALL_b37_QC_grm_chr$chr --thread-num 5 #--keep
	done
	##make the input file that lists all the grm
	for chr in $(seq 1 22); do 
	echo "$workdir/gcta/clean_ALL_b37_QC_grm_chr"$chr ; done > $workdir/gcta/cleanmulti_ALL_grm.txt
	# adjusting the grm
	$gcta --mgrm $workdir/gcta/cleanmulti_ALL_grm.txt --grm-adj 0 --grm-cutoff 0.025 --make-grm --out $workdir/gcta/clean_b37_QC_grm_ALL
	#Then calculating the PC's
	$gcta --grm $workdir/gcta/clean_b37_QC_grm_ALL  --pca 20 --out $workdir/gcta/clean_b37_QC_PC_ALL
	#Then running the actual heritability estimate with prevalence equal to the reported Disease prevalence for the country. 
	$gcta --reml --grm $workdir/gcta/clean_b37_QC_grm_ALL --pheno $workdir/gcta/clean_phenoforgcta_ALL.txt --qcovar $workdir/gcta/clean_b37_QC_PC_ALL.eigenvec --prevalence $ALL_prevalence --out $workdir/gcta/clean_b37_QC_actual_ALL

"


#使用ANNOvar进行基因注释
#使用plink将ped文件转成vcf文件
mkdir annovar
#plink --bfile genotype/clean.sex.hwe.rmrelated.rmoutlier --recode vcf --out annovar/clean
#将vcf格式转化为annovar intput格式
$annovar/convert2annovar.pl -format vcf4 -allsample -withfreq  annovar/clean.vcf > clean.avinput











