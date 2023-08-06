本仓库曾经使用过现已经废弃，usci 为最新分支，后续会在近期可能发表的文章补充材料中，重新整理相关功能及计算结果到新的仓库。（所有功能均已重新架构融入http://www.atcgorder.top:8082/v1/user/login分析平台id:atcgorder,pw:admin专利正在申请中）
本仓库的大量函数模块仍然可以使用：
NGS\Analysis下
usedadiPy2_7\GenerateSNPfilefromvcftable.py
usedadiPy2_7\bootstrapdadisimulation.py
为实施群体历史dadi模拟代码，GenerateSNPfilefromvcftable.py（GenerateSNPfilefromvcffile.py）从提供的vcf文件中按照规则提取snp并生成dadi输入文件。bootstrapdadisimulation.py按照规则随机抽取进行重复模拟，调用dadicode.py、dadicode_split.py等并执行指定模型。实际调用dadi软件是修改过源代码的版本以满足将图像数据提取出来，得到综合统计结果。
CalculateKaKsusingPaml.py使用见‘利用paml批量计算基因的dnds及记录.pdf’
one2oneOrtholog.py为趋同进化分析'趋同进化.docx'
ancient seletion早期选择分析.docx 实现了‘A Draft Sequence of the Neandertal Genome’sience2010中Richard E. Green开发的探测分化发生时,sister species当中某一支相对另一支受到选择的基因组片段。
全基因组滑窗计算统计量dxy,df,fst,hp,pi等，使用Detectsignalacrossgenome_master.py（计算逻辑调用Calculators.py）

......



版权所有，翻版必究！
