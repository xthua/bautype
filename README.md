Author:     Xiaoting Hua at 20190806

Email:      XiaotingHua@zju.edu.cn

institute:  Key laboratory of Microbiol technology and Bioinformatics of Zhejiang Province

This program is the  pipeline of LPS or CPS typing for Acinetobacter baumannii


Run:

    bautype is a python3.X script, compatible with Linux and windows.
    1.BLAST is required, it shoud be added in environment variable.
      read requirements.txt , check your dependecies or just run: pip install -r requirements.txt
    2.Currently, it only process FASTA format. You can just run like following:
      for K type: python AB_LPSandCPS_typeing.py  -o outdir -m tmpdir -r reference_data/Abaumannii_KL_reference -i fastaFile
                  python AB_LPSandCPS_typeing.py  -o outdir -m tmpdir -r reference_data/Abaumannii_KL_reference -d indir
      for O type: python AB_LPSandCPS_typeing.py  -o outdir -m tmpdir -r reference_data/Abaumannii_OCL_reference -i fastaFile
                  python AB_LPSandCPS_typeing.py  -o outdir -m tmpdir -r reference_data/Abaumannii_OCL_reference -d indir
      
Output:

    output_match_result.txt        LPS or CPS typing result for your input
    locus_gene.coverage.svg/png    graph for the gene coverage
