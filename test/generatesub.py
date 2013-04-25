import io
f = io.open('rand.sub','wb')

f.write('getenv = TRUE\n')
f.write('Universe = vanilla\n')
f.write('Executable = randtest.py\n')

for i in ['b','q','w']:
    for k in ['b','w','z']:
        for r in range(1,51):
            fname = i+k+str(r)
            f.write('arguments = '+fname+'\n')
            f.write('Output = randresults/outs/'+fname+'.out\n')
            f.write('Error = randresults/errs/'+fname+'.err\n')
            f.write('Log = randresults/logs/'+fname+'.log\n')
            f.write('queue\n\n')

f.close()
