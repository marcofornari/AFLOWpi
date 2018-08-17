import get_K_grid_fft
import numpy as np


b_vectors = np.array([[1,0,0,],
                      [0,1,0,],
                      [0,0,1,]])

nk3 = 128
nk2 = 128
nk1 = 128

grid,ki,nkt,idk = get_K_grid_fft.get_K_grid_fft(nk1,nk2,nk3,b_vectors)

grid=grid


grid = np.reshape(grid,(3,nk1,nk2,nk3),order='C')

grid =np.transpose(grid,axes=(1,2,3,0))%1.0

np.set_printoptions(precision=3,suppress=True,formatter={'all':lambda x: '% 3.3f'%x})

print grid.shape

#grid = np.roll(grid,grid.shape[0]/2,axis=0)
#grid = np.roll(grid,grid.shape[1]/2,axis=1)
#grid = np.roll(grid,grid.shape[2]/2,axis=2)

# for i in xrange(grid.shape[0]):
#     print ''
#     for j in xrange(grid.shape[1]):
#         print grid[i,j,0,:],' ',grid[i,j,1,:],' ',grid[i,j,2,:],' ',grid[i,j,3,:]#,' ',grid[i,j,4,:],' ',grid[i,j,5,:]

print ''
print ''


new_grid = np.zeros((nk1+1,nk2+1,nk3+1,3))
new_grid[:-1,:-1,:-1]=grid

#new_grid[-1,-1,-1] = new_grid[0,0,0]


new_grid[-1, :, :] = new_grid[0,:,:]
new_grid[ :,-1, :] = new_grid[:,0,:]
new_grid[ :, :,-1] = new_grid[:,:,0]

# new_grid[ :,-1,-1] = new_grid[:,0,0]
# new_grid[-1, : -1] = new_grid[0,:,0]
# new_grid[-1,-1, :] = new_grid[0,0,:]




# new_grid[-1, :, :] = new_grid[0,:,:]
# new_grid[ :,-1, :] = new_grid[:,0,:]
# new_grid[ :, :,-1] = new_grid[:,:,0]

new_grid%=1.0






eps=0.0
a=np.linspace(0.0-eps,1.0+eps,nk1+1,endpoint=True)
b=np.linspace(0.0-eps,1.0+eps,nk2+1,endpoint=True)
c=np.linspace(0.0-eps,1.0+eps,nk3+1,endpoint=True)
X,Y,Z =np.meshgrid(a,b,c,indexing='ij')

test_arr=np.array([X,Y,Z])

test_arr =np.transpose(test_arr,axes=(1,2,3,0))%1.0



print np.all(test_arr[:,:,:]==new_grid[:,:,:])



# for i in xrange(new_grid.shape[0]):
#     print ''
#     for j in xrange(new_grid.shape[1]):
#         print new_grid[i,j,0,:],' ',new_grid[i,j,1,:],' ',new_grid[i,j,2,:],' ',new_grid[i,j,3,:]#,' ',new_grid[i,j,4,:],' ',new_grid[i,j,5,:]
# print ''
# print ''
# for i in xrange(grid.shape[0]):
#     print ''
#     for j in xrange(grid.shape[1]):
#         print grid[i,j,0,2],' ',grid[i,j,1,2],' ',grid[i,j,2,2],' ',grid[i,j,3,2]

print ''
print ''



eps=0.0
a=np.linspace(0.0-eps,1.0+eps,nk1+1,endpoint=True)
b=np.linspace(0.0-eps,1.0+eps,nk2+1,endpoint=True)
c=np.linspace(0.0-eps,1.0+eps,nk3+1,endpoint=True)
X,Y,Z =np.meshgrid(a,b,c,indexing='ij')

test_arr=np.array([X,Y,Z])

test_arr =np.transpose(test_arr,axes=(1,2,3,0))%1.0

print new_grid.shape

# print np.all(test_arr[:,:,:,1]==new_grid[:,:,:,1])
# for i in xrange(test_arr.shape[0]):
#     print ''
#     for j in xrange(test_arr.shape[1]):
#         print test_arr[i,j,0,:],' ',test_arr[i,j,1,:],' ',test_arr[i,j,2,:],' ',test_arr[i,j,3,:]#,' ',test_arr[i,j,4,:]
