################################################################################
### Package imports ###
using DifferentialEquations
using Plots
GR()

################################################################################
### Functions ###
function Motion!(d,PX,PY,PZ,VX,VY,VZ,G,N,M,T)
    A = zeros(N)
    AN = N*3
    for n = 1:N
        for i = 1:3 # spatial dimension
            d[3*(n-1)+i] = d[AN+3*(n-1)+i] # velocity of object n in dimension i
        end
        d[AN+3*(n-1)+1] = sum((G*M[ii]*(PX[ii]-PX[1]))/((PX[ii]-PX[1])^2+(PY[ii]-PY[1])^2+(PX[ii]-PX[1])^2)^(1.5)), ii = 1:3 if ii != 1) # accelration
        d[AN+3*(n-1)+2] = sum((G*M[ii]*(PX[ii]-PY[2]))/((PX[ii]-PX[1])^2+(PY[ii]-PY[1])^2+(PX[ii]-PX[1])^2)^(1.5)), ii = 1:3 if ii != 1) # accelration
    end

    # Motion!(du,u,p,t)
    r,w = u
    α,β,γ,δ=p
    du[1]= dr = α*r - β*r*w
    du[2]= dw = γ*r*w - δ*w
end






################################################################################
### Parameters ###
G = 1 # gravitation constant, normalized

N = 3 # number of objects

M = [1,0.75,1.25] # Object masses

T = (0.0,100.0) # time horizon

PX0 = [10.0,0.0,0.0] # Position
PY0 = [0.0,-10.0,0.0]
PZ0 = [0.0,0.0,10.0]

VX0 = [1.0,-1.0,0.0] # Velocity
VY0 = [0.5,1.0,-0.3]
VZ0 = [0.3,-0.3,0.5]
