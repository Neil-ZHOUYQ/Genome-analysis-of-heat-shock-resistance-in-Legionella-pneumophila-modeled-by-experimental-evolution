from Bio import Entrez, SeqIO

# 邮箱地址
Entrez.email = "lwy020728@gmail.com"  

# 想要下载的参考序列号码
accession_number = str(input("Please enter the ref_number: "))
# "NC_002516.2"

try:
    # 获取数据
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    record = handle.read()

    # 将数据保存为本地文件
    with open(f"{accession_number}.fasta", "w") as output_file:
        output_file.write(record)

    print(f"下载完成：{accession_number}.fasta")

except Exception as e:
    print(f"发生错误：{e}")
