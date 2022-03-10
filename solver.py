import numpy as np
import shutil
import os

def pressure_poisson(p, dx, dy, rho, dt, u, v, pBCs):
    pn = np.empty_like(p)
    pn = p.copy()
    nit = 50   #pseudo-time steps in each timestep
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + 
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                          (2 * (dx**2 + dy**2)) -
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * ((rho * (1 / dt * 
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / 
                     (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                      2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                           (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                          ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))))


        pLeft=pBCs[0]
        pRight=pBCs[1]
        pTop=pBCs[2]
        pBottom=pBCs[3]

        pressureBC(pLeft,pRight,pTop,pBottom,p,dx,dy)
        #p[:, -1] = 0 #outlet pressure of 0
        #p[0, :] = p[1, :]   # dp/dy = 0 at bottom
        #p[:, 0] = 1 # inlet pressure of 1
        #p[-1, :] = p[-2,:]  # dp/dx = 0 at top
        
    return p



def flow_solver(nt, u, v, dt, ds, dx, dy, p, rho, nu, pBCs, uBCs, vBCs):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    #b = np.zeros((ny, nx))
    i=0
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        i=i+1
        #b = bf(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, rho, dt, u, v, pBCs)
        
        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        
        velmag=(u**2+v**2)**0.5

        if ((n+1)%ds==0): #saving the results into .csv files
            np.savetxt('u'+str(i)+'.csv', u, delimiter=',')
            np.savetxt('v'+str(i)+'.csv', v, delimiter=',')
            np.savetxt('p'+str(i)+'.csv', p, delimiter=',')
            np.savetxt('Umag'+str(i)+'.csv', velmag, delimiter=',')
            

        #boundary conditions for u

        uLeft=uBCs[0]
        uRight=uBCs[1]
        uTop=uBCs[2]
        uBottom=uBCs[3]
        XVelBC(uLeft,uRight,uTop,uBottom,u,dx,dy)

        #u[0, :]  = 0    #bottom
        #u[-1, :] = 0    #top
        #u[:, 0]  = u[:, 1]   #left
        #u[:, -1] = u[:,-2]    #right
        #boundary conditions for v
        vLeft=vBCs[0]
        vRight=vBCs[1]
        vTop=vBCs[2]
        vBottom=vBCs[3]
        YVelBC(vLeft,vRight,vTop,vBottom,v,dx,dy)
        #v[0, :]  = 0    #same as u
        #v[-1, :] = 0
        #v[:, 0]  = 0
        #v[:, -1] = 0
                
    return u, v, p

def pressureBC(a,b,c,d,p,dx,dy):
    #left
    if a[0]=='N':
        p[:, 0]=a[1]*dx+p[:, 1]
    elif a[0]=='D':
        p[:,0]=a[1]
    else:
        return('Please enter a valid BC type')

    #right
    if b[0]=='N':
        p[:, -1]=b[1]*dx+p[:, -2]
    elif b[0]=='D':
        p[:, -1]=b[1]
    else:
        return('Please enter a valid BC type')

    #top
    if c[0]=='N':
        p[-1,:]=c[1]*dy+p[-2, :]
    elif c[0]=='D':
        p[-1,:]=c[1]
    else:
        return('Please enter a valid BC type')

    #bottom
    if d[0]=='N':
        p[0,:]=d[1]*dy+p[1, :]
    elif d[0]=='D':
        p[0,:]=d[1]
    else:
        return('Please enter a valid BC type')

def XVelBC(a,b,c,d,u,dx,dy):
    #left
    if a[0]=='N':
        u[:, 0]=a[1]*dx+u[:, 1]
    elif a[0]=='D':
        u[:,0]=a[1]
    else:
        return('Please enter a valid BC type')

    #right
    if b[0]=='N':
        u[:, -1]=b[1]*dx+u[:, -2]
    elif b[0]=='D':
        u[:, -1]=b[1]
    else:
        return('Please enter a valid BC type')

    #top
    if c[0]=='N':
        u[-1,:]=c[1]*dy+u[-2, :]
    elif c[0]=='D':
        u[-1,:]=c[1]
    else:
        return('Please enter a valid BC type')

    #bottom
    if d[0]=='N':
        u[0,:]=d[1]*dy+u[1, :]
    elif d[0]=='D':
        u[0,:]=d[1]
    else:
        return('Please enter a valid BC type')

def YVelBC(a,b,c,d,v,dx,dy):
    #left
    if a[0]=='N':
        v[:, 0]=a[1]*dx+v[:, 1]
    elif a[0]=='D':
        v[:,0]=a[1]
    else:
        return('Please enter a valid BC type')

    #right
    if b[0]=='N':
        v[:, -1]=b[1]*dx+v[:, -2]
    elif b[0]=='D':
        v[:, -1]=b[1]
    else:
        return('Please enter a valid BC type')

    #top
    if c[0]=='N':
        v[-1,:]=c[1]*dy+v[-2, :]
    elif c[0]=='D':
        v[-1,:]=c[1]
    else:
        return('Please enter a valid BC type')

    #bottom
    if d[0]=='N':
        v[0,:]=d[1]*dy+v[1, :]
    elif d[0]=='D':
        v[0,:]=d[1]
    else:
        return('Please enter a valid BC type')


def clearResults():
    fileList=os.listdir('./')
    
    if os.path.isdir('./Results'):
        shutil.rmtree('./Results')
        os.mkdir('Results')
    else:
        os.mkdir('Results')
    
    for file in fileList:
        if file.endswith('.csv'):
            shutil.move(file,'Results/.')

