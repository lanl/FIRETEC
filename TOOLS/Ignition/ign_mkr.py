#write new ignite.dat files
import numpy as np
import math
import os

def main():
    name_dict = {0:'No fire',1:'Line fire',
                2:'Dot fire (aerial)',
                3:'Multiple line fire (atv)',
                22: 'DEBUGGING CHOICE - JULIA'}
    print('Hit ctrl+c to end program at any time......')
    print('Ignition Patterns: ')
    print(name_dict.keys())
    for key in list(name_dict.keys())[0:4]:
        print(key,' = ',name_dict[key])

    igntype = tryResponse('Enter ignition pattern choice: ')

    while True:
        try:
            print(name_dict[igntype]+' selected')
        except KeyError:
            print('Please select a valid ignition strategy: ')
            print('Ignition Patterns: ')
            for key in list(name_dict.keys())[0:4]:
                print(key,' = ',name_dict[key])
            igntype = tryResponse('Enter ignition pattern choice: ')
            continue
        else:
            break

    if igntype>0:
        rr = tryResponse('Enter ramp rate (K/s), ex=350 : ')
        tt = tryResponse('Enter target temperature (K), ex=1000 : ')
        strt_time = tryResponse('Enter start time (s) : ')
        if igntype==1: #LINE FIRE INPUTS
            minx      = tryResponse('Enter min x location (cells) : ')
            maxx      = tryResponse('Enter max x location (cells) : ')
            miny      = tryResponse('Enter min y location (cells) : ')
            maxy      = tryResponse('Enter max y location (cells) : ')
            ign_all   = lineIgn(minx,maxx,miny,maxy,strt_time)
            numignpts = len(ign_all[:,1])
            writeIgn(igntype,numignpts,tt,rr,ign_all,name='ignite.dat')
        if igntype==2: #AERIAL INPUTS
            staggered  = tryResponse('Enter staggered (0=NO, 1=YES) : ')
            ltr        = tryResponse('Enter drop order left-to-right (0=NO, 1=YES) : ')
            print('Ignition drop size currently only produces squares! ')
            print('Dropsize is length or width of single ignition point')
            dropsize   = tryResponse('SEE ABOVE -- Enter ignition drop size (cells) : ') 

            nrcolumns  = tryResponse('Enter number of ignition rows : ')
            nrrows     = tryResponse('Enter number of ignition drops in a row : ')
            startcellx = tryResponse('Enter start x location (cells) : ')
            startcelly = tryResponse('Enter start y location (cells) : ')
            spacingx   = tryResponse('Enter spacing of ignition drops in x (cells) : ')
            spacingy   = tryResponse('Enter spacing of ignition drops in y (cells) : ')
            turntime   = tryResponse('Enter aerial turn time (s), ex=10 : ')
            speed      = tryResponse('Enter speed of drops (m/s), ex=15.6 : ') #aerial craft speed m/s; drone; 35mph
            di         = tryResponse('Enter grid spacing (dx/dy), ex=2m :')

            ign_all = aerialIgn(staggered,ltr,dropsize,nrrows,nrcolumns,startcellx,startcelly,
                                spacingx,spacingy,strt_time,turntime,speed,di)
            numignpts = len(ign_all[:,1])
            writeIgn(igntype,numignpts,tt,rr,ign_all,name='ignite.dat')

        if igntype==3: #ATV INPUTS
            dx              = tryResponse('Enter horizontal cell resolution, ex: 2m : ')
            atvspeed        = tryResponse('Enter speed of ATV (m/s), ex: 7.5mph = 3.353m/s : ')
            ltr             = tryResponse('Enter driving order left-to-right (0=NO, 1=YES) : ')
            ncols           = tryResponse('Enter total number of ATV lines : ')
            
            #natv         = tryResponse('Enter set number of atvs to run at once, ex: 3 : ')  [ADD ABILITY TO DO MULTIPLE LINES IN A GROUP]
            #staggered = 4#int(input('Enter stagger distance (cells) between lines, ex: 2 : ')) [ ADD STAGGER BETWEEN MULTIPLE LINES]
            
            flip = tryResponse('Do you want to flip the direction of ATV lines after each successive set of ATVs (0 = NO, 1 = YES)? " ')
            
            dash = tryResponse('Enter the pattern you want to make (0 = Line, 1 = Dash, 2 = Alternating dash) : ')
            if dash > 0:
                dash_sz = tryResponse('Enter the dash size (cells) : ')
                space_sz = tryResponse('Enter the space between dashes (cells) : ')
            else:
                dash_sz = 0
                space_sz = 0
                
            startcellx = tryResponse('Enter start x location of first line (cells) : ')
            startcelly = tryResponse('Enter start y location of first line (cells) : ')

            length    = tryResponse('Enter the maximum length (cells) of the ATV line  : ')
            spacing   = tryResponse('Enter spacing (cells) between the ATV lines : ')
            turntime  = tryResponse('Enter ATV turn time (s) : ')

            ign_all = atvIgn(dx, atvspeed, ltr, ncols, startcellx, startcelly, 
                             length, spacing, strt_time, turntime,dash,flip,
                             dash_sz,space_sz)
            numignpts = len(ign_all[:,1])
            writeIgn(igntype,numignpts,tt,rr,ign_all)

        if igntype==22: #JULIA'S QUICWRITE DEBUGGER DON'T TOUCH
            #QUICK WRITE FOR TESTING - use on 100 x 100 cell domain
            #NO IGN
            writeIgn(0,0,0,0,np.zeros((2,2)),name='ignite_none.dat')
            #LINE FIRE
            ign_all   = lineIgn(25,30,25,30,strt_time)
            numignpts = len(ign_all[:,1])
            writeIgn(1,numignpts,tt,rr,ign_all,name='ignite_line.dat')
            #AERIAL FIRE
            ign_all = aerialIgn(1,1,3,10,3,25,25,
                                3,3,strt_time,5,15.6,2)
            numignpts = len(ign_all[:,1])
            writeIgn(2,numignpts,tt,rr,ign_all,name='ignite_aerial.dat')
            #ATV FIRE
            ign_all = atvIgn(2, 3.353, 1, 3, 25, 25, 
                             50, 3, strt_time, 10,0,1,
                             0,0)
            numignpts = len(ign_all[:,1])
            writeIgn(3,numignpts,tt,rr,ign_all,name='ignite_atv.dat')
            
    else: 
        numignpts = 0
        tt        = 0
        rr        = 0
        writeIgn(igntype,numignpts,tt,rr,np.zeros((2,2)),name='ignite.dat')
    return 0


####### functions ########
def lineIgn(minx,maxx,miny,maxy,strt_time):
    minx = int(minx)
    miny = int(miny)
    maxx = int(maxx)
    maxy = int(maxy)
    ign_all = np.zeros((((maxx-minx)+1)*((maxy-miny)+1),3))
    k = 0
    for i in range( (maxx-minx)+1 ):
        for j in range( (maxy-miny)+1 ):
            ign_all[k,0] = i+minx
            ign_all[k,1] = j+miny
            ign_all[k,2] = np.round(strt_time,2)
            k+=1
    return ign_all

def aerialIgn(staggered,ltr,dropsize,nrrows,nrcolumns,startcellx,startcelly,
              spacingx,spacingy,starttime,turntime,speed,di):

    nrcolumns  = int(nrcolumns)
    nrrows     = int(nrrows)
    dropsize   = int(dropsize)
    spacingx   = int(spacingx+dropsize) #cells
    spacingy   = int(spacingy+dropsize) #cells

    #=============================================================
    if staggered:
        if isinstance(nrcolumns/2, int):
            ignpoints = np.zeros(((dropsize**2 * (nrrows*nrcolumns/2)) + (dropsize**2 * ((nrrows-1)*nrcolumns/2)),2), dtype=int)
        else:
            staggcols = int(np.floor(nrcolumns/2))
            ignpoints = np.zeros(((dropsize**2 * (nrrows*(nrcolumns-staggcols))) + (dropsize**2 * ((nrrows-1)*staggcols)),2), dtype=int)
    if not staggered:
        ignpoints = np.zeros((dropsize**2 * (nrrows*nrcolumns),2), dtype=int)

    igntime = np.zeros( len(ignpoints) , dtype=np.float32)

    k   = 0
    b   = 1
    fpc = 0 #first pass column
    fpr = 0 #first pass row

    if ltr:
        for columncount in range(nrcolumns):
            for rowcount in range(nrrows):
                if staggered and (columncount%2==1) and (rowcount==nrrows-1):
                    break
                for bb in range(0,dropsize):
                    x = int(startcellx+bb+(spacingx*columncount))
                    for bbb in range(0,dropsize):
                        if staggered and (columncount%2==1):
                            y = int(startcelly+bbb+(spacingy*rowcount)+int(spacingy/2))
                        else:
                            y = int(startcelly+bbb+(spacingy*rowcount))
                        ignpoints[k,0] = x
                        ignpoints[k,1] = y
                        igntime[k] = np.round(starttime+(fpr*spacingy*di/speed)+(fpc*turntime),2)
                        k+=1 #iterate through all ignpoints
                b+=1 #keep track of balls dropped
                fpr+=1
            fpc+=1
    else:    
        for columncount in range(nrcolumns):
            for rowcount in range(nrrows):
                if staggered and (columncount%2==1) and (rowcount==nrrows-1):
                    break
                for bb in range(0,dropsize):
                    x = int(startcellx-bb-(spacingx*columncount))
                    for bbb in range(0,dropsize):
                        if staggered and (columncount%2==1):
                            y = int(startcelly-bbb-(spacingy*rowcount)-int(spacingy/2))
                        else:
                            y = int(startcelly-bbb-(spacingy*rowcount))
                        ignpoints[k,0] = x
                        ignpoints[k,1] = y
                        igntime[k] = np.round(starttime+(fpr*spacingy*di/speed)+(fpc*turntime),2)
                        k+=1 #iterate through all ignpoints
                b+=1 #keep track of balls dropped
                fpr+=1
            fpc+=1

    stg = 0
    cnt = 0
    switch = 0
    if staggered:
        ign_step = np.zeros(nrcolumns).astype(int)
        for kk in range(1, nrcolumns):
            if ((kk%2)==0):
                switch = dropsize**2
            ign_step[kk] = ign_step[kk-1]+nrrows*dropsize**2-switch
            switch = 0
    else:
        ign_step = np.arange(0,len(ignpoints[:,0]),nrrows*dropsize**2)

    for w in ign_step:
        if staggered:
            stg=1
        if ((cnt%2)==1): 
            arr = igntime[w:w+(nrrows-stg)*dropsize**2]
            arr = arr[::-1]
            igntime[w:w+(nrrows-stg)*dropsize**2] = arr
        cnt +=1

    ign_all = np.zeros((len(igntime), 3)).astype(object)
    ign_all[:,0] = ignpoints[:,0].astype(int)
    ign_all[:,1] = ignpoints[:,1].astype(int)
    ign_all[:,2] = np.round(igntime,2)

    for w in ign_step:
        if staggered:
            stg=1
        if ((cnt%2)==1): 
            arr = ign_all[w:w+(nrrows-stg)*dropsize**2,:]
            arr = arr[::-1]
            ign_all[w:w+(nrrows-stg)*dropsize**2,:] = arr
        cnt +=1
    
    ign_all = ign_all[np.logical_not(np.logical_or( ign_all[:,0] < 0, ign_all[:,1] < 0 )) ]
    return ign_all

#######
def atvIgn(dx, atvspeed, ltr, ncols, startcellx, startcelly, length, spacing, strt_time, turntime,dash,flip,dash_sz,space_sz):
    #lines

    time = strt_time
    ncols = int(ncols)
    length = int(length)
    dash_sz = int(dash_sz)
    xl = []
    yl = []
    tm = []
    if ltr: 
        if dash == 0: #Line patterns
            for col in range(ncols):
                xloc = startcellx + col * (spacing + 1)
                
                for yuh in range(length):
                    if flip == 0: #dont flip direction
                        yloc = startcelly + yuh
                    if flip == 1: #flip direction
                        if col % 2 == 0:
                            yloc = startcelly + yuh
                        else:
                            yloc = startcelly + length - (1 + yuh)
                            
                    #time
                    if yuh == 0 and col == 0:
                        time = time
                    elif yuh == 0 and col > 0:
                        time += turntime
                    else:
                        time += (dx / atvspeed)

                    xl.append(xloc) ; yl.append(yloc) ; tm.append(np.round(time,2))

        dash = int(dash)
        if dash == 1: #Dash patterns
            ndash = int(math.floor((length) / (dash_sz + space_sz))) #find number of dashes that will fit

            for col in range(ncols):
                xloc = startcellx + col * (spacing + 1 )
                
                for da in range(ndash):
                    for spa in range(dash_sz):
                        if flip == 0:
                            yloc = startcelly + spa + da * (space_sz + dash_sz)
                        if flip == 1:
                            if col % 2 == 0: #if even
                                yloc = startcelly + spa + da * (space_sz + dash_sz)
                            else: #if odd
                                yloc = ndash * (space_sz + dash_sz) - space_sz - (spa + da * (space_sz + dash_sz))

                        #timing
                        if col == 0 and da == 0 and spa == 0:
                            time = time
                        elif col > 0 and da == 0 and spa == 0:
                            time += turntime
                        elif spa == 0:
                            time += (space_sz * dx / atvspeed)
                        else:
                            time += (dx / atvspeed)
                            
                        xl.append(xloc) ; yl.append(yloc) ; tm.append(np.round(time,2))
        
    else:
        if dash == 0: #Line patterns
            for col in range(ncols):
                xloc = startcellx - col * (spacing + 1)
                
                for yuh in range(length):
                    if flip == 0: #dont flip direction
                        yloc = startcelly + yuh
                    if flip == 1: #flip direction
                        if col % 2 == 0:
                            yloc = startcelly + yuh
                        else:
                            yloc = startcelly + length - (1 + yuh)
                            
                    #time
                    if yuh == 0 and col == 0:
                        time = time
                    elif yuh == 0 and col > 0:
                        time += turntime
                    else:
                        time += (dx / atvspeed)

                    xl.append(xloc) ; yl.append(yloc) ; tm.append(np.round(time,2))

        dash = int(dash)
        if dash == 1: #Dash patterns
            ndash = int(math.floor((length) / (dash_sz + space_sz))) #find number of dashes that will fit

            for col in range(ncols):
                xloc = startcellx - col * (spacing + 1 )
                
                for da in range(ndash):
                    for spa in range(dash_sz):
                        if flip == 0:
                            yloc = startcelly + spa + da * (space_sz + dash_sz)
                        if flip == 1:
                            if col % 2 == 0: #if even
                                yloc = startcelly + spa + da * (space_sz + dash_sz)
                            else: #if odd
                                yloc = ndash * (space_sz + dash_sz) - space_sz - (spa + da * (space_sz + dash_sz))

                        #timing
                        if col == 0 and da == 0 and spa == 0:
                            time = time
                        elif col > 0 and da == 0 and spa == 0:
                            time += turntime
                        elif spa == 0:
                            time += (space_sz * dx / atvspeed)
                        else:
                            time += (dx / atvspeed)
                            
                        xl.append(xloc) ; yl.append(yloc) ; tm.append(np.round(time,2))

    ign_all = np.zeros((len(xl),3))
    ign_all[:,0] = xl
    ign_all[:,1] = yl
    ign_all[:,2] = tm
    ign_all = ign_all[np.logical_not(np.logical_or( ign_all[:,0] < 0, ign_all[:,1] < 0 )) ]
    return ign_all

def writeIgn(igntype,numignpts,tt,rr,ign_all,name):
    pf = os.getcwd()
    fname = pf+'/'+name
    if os.path.isfile(fname):
        os.remove(fname)
    f = open(fname,'w')
    print('Writing ignition file')

    #WRITE HEADERS!
    f.write('&ignitelist\n')
    if igntype==0:
        f.write('  nIgn=0\n')
        f.write('  rampRate=0\n')
        f.write('  targetTemp=300\n')
        f.write('/ \n')
        f.close()
    else:
        f.write('  nIgn='      +str(int(numignpts))+'\n')
        f.write('  rampRate='  +str(int(rr))       +'\n')
        f.write('  targetTemp='+str(int(tt))       +'\n')
        f.write('/ \n')
        #format: x cell, y cell, time
        for ii in range(numignpts):
            f.write('    '+str(int(ign_all[ii,0]))+'    '+str(int(ign_all[ii,1]))+'    '+str(1)+'    '+str(np.round(ign_all[ii,2],2))+'\n')
        f.close()
    return

def tryResponse(input_str):

    while True:
        try:
            val = float(input(input_str))
        except ValueError:
            print("Error, please enter numerical value")
            continue
        else:
            #successfully parsed!
            break
    return val

main()
