if __name__ == '__main__':
    filename = r'E:\01_S915\02_PLA2R\00_SvnRepo\02_PLA2R表达质粒信息与验证\01_质粒序列\01_NCBI序列信息\PLA2R1-mRNA.fasta'
    with open(filename) as f:
        content = f.readlines()
        gene_info = content[0]
        sequence = ''.join(content[1:]).replace('\n','')

    from DNA import DNA

    s = DNA(sequence)
    print(s.translate())
    print(s.orf)
        
