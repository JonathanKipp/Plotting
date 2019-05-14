from datetime import date
today = date.isoformat(date.today())
prefix = "/Users/kipp/STB/bashfiles/"
seedfname = "kipp_batch"
files = 4
lambdas = 6 
tsos = 6
#filepath = "/p/home/jusers/kipp1/jureca/STB/inis/lambda_tso/path_rel_G-K-Kprime_anticol_theta_scananticol"
#filepath2 = "/p/home/jusers/kipp1/jureca/STB/inis/lambda_tso/path_rel_G-K-Kprime_anticol_theta_scananticolflip"
filepath = "/p/home/jusers/kipp1/jureca/STB/inis/path_rel_G-K-Kprime_anticol_theta_scananticol"
filepath2 = "/p/home/jusers/kipp1/jureca/STB/inis/path_rel_G-K-Kprime_anticol_theta_scananticolflip"
#keydict = {account:"jiff40",nodes:1,ntasks:68,ntasks-per-node:68,output:"mpi-out.%j",error:"mpi-err.%j",time:"04:00:00",partition:"develbooster",mail-user:"j.kipp@fz-juelich.de",mail-type:"BEGIN,END,FAIL"}
for j in range(files):
    for k in range(lambdas):
        for l in range(tsos):
            with open(prefix + seedfname,'rt') as fin:
                with open(prefix + today + "/" + seedfname + "_{:02}_lambda_{:02}_t_so_{:02}".format(j,k,l), "wt") as fout:	
                    for line in fin:
                        fout.write(line)
                    writestring = "srun /p/home/jusers/kipp1/jureca/STB/stb.x " + filepath +  "_{:02}_lambda_{:02}_t_so_{:02}.cfg\n".format(j,k,l)
                    writestring2 = "srun /p/home/jusers/kipp1/jureca/STB/stb.x " + filepath2 +  "_{:02}_lambda_{:02}_t_so_{:02}.cfg".format(j,k,l)
                    fout.write(writestring)
                    fout.write(writestring2)
                fout.close()