'''
Test randomly problem with positive definite Hessian
'''
#from subprocess import call
import os
def sprandsym(n,den,rc,mtype = 1,dataname = 'tmp'):
#	curdir = os.path.dirname(os.path.realpath(__file__)).replace('random.pyc','')
	curdir = os.path.abspath(__file__).replace('random.pyc','')
	curdir = curdir.replace('random.py','')
	print curdir
	st = 'matlab -nodisplay -r '+'"addpath '+curdir+';random_bqp('+ str(n) +','+ str(den)+','+str(rc)+','+str(mtype)+','+"'" +str(dataname)+"'"+');exit;"'
	print st
	os.system(st)

# Call as function
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        mtype = 1
    else:
        mtype = sys.argv[4]
    if len(sys.argv) < 6:
        dname = 'tmp'
    else:
        dname = sys.argv[5]

    print len(sys.argv)
    n = sys.argv[1]
    print n
    den = sys.argv[2]
    print den
    rc = sys.argv[3]
    print rc
    print mtype
    print dname
    sprandsym(n,den,rc,mtype,str(dname))
