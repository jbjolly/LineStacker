import numpy as np
import sys

def ProgressBar(n, total):
    """
        show progress of a process
        Parameters
        ---------
        n
            current step
        total
            total number of steps
    """
    if total>100:
        g=int(10*n/total)
        if n in np.arange(0,total, int(total/100.)) and n!=total-1:
            #done = ("\033[1;36m"+'#'+'\x1b[0m')*(g)
            #todo = '-'*(10-g)
            done=("\033[0;36;46m"+" "+'\x1b[0m')*g
            todo = ' '*(10-g)
            s = '|{0}|'.format(done+todo)
            s+="\033[1;33m"+' LOADING '+'\x1b[0m'
            s+="\033[0;32m"+str(int((float(n)/total)*100))+'% '+'\x1b[0m'
            s = '\r'+s
            print s,
            sys.stdout.flush()
            #return
        if n==total-1:
            g=10
            done=("\033[0;36;46m"+" "+'\x1b[0m')*g
            #todo = '-'*(10-g)
            s = '|{0}|'.format(done)
            s+="\033[1;33m"+' LOADING '+'\x1b[0m'
            s+='\033[1;31m'+ ' 100%' +'\x1b[0m'+ '\n\n'
            s = '\r'+s
            print s,
            sys.stdout.flush()
    else:
        if n==total-1:
            g=10
            done=("\033[0;36;46m"+" "+'\x1b[0m')*g
            #todo = '-'*(10-g)
            s = '|{0}|'.format(done)
            s+="\033[1;33m"+' LOADING '+'\x1b[0m'
            s+='\033[1;31m'+ ' 100%' +'\x1b[0m'+ '\n\n'
            s = '\r'+s
            print s,
            sys.stdout.flush()
        else:
            ratio=float(n)/total
            g=int(ratio*10)
            done=("\033[0;36;46m"+" "+'\x1b[0m')*g
            todo = ' '*(10-g)
            s = '|{0}|'.format(done+todo)
            s+="\033[1;33m"+' LOADING '+'\x1b[0m'
            s+="\033[0;32m"+str(int((float(n)/total)*100))+'% '+'\x1b[0m'
            s = '\r'+s
            print s,
            sys.stdout.flush()
