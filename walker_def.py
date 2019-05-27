import numpy as np

def walker(scale,N):
    '''
    INPUTS:
    scale: scaling of the quantity
    N: number of data points
    '''
    
    # probability to move up or down 
    prob = [0.05, 0.95] 

    # step scale factor
    step_scale = 0.005
    
    # starting position 
    start = 0
    positions = [start] 
    
    # creating the random points 
    rr = np.random.random(N) 
    downp = rr < prob[0] 
    upp = rr > prob[1] 
    
    for idownp, iupp in zip(downp, upp): 
        down = idownp and positions[-1] > -1
        up = iupp and positions[-1] < 1
        positions.append(positions[-1] - step_scale*down + step_scale*up) 
        
    positions_scaled = np.multiply(scale,positions)
        
    return positions_scaled
