################################################################################
### JULIA PACKAGE IMPORTS ###
using Plots
using Colors
using DelimitedFiles

################################################################################
### CONTROL FLOW ###

# Source data, if any, also the save location for the next run
Xfilename = "LastX.txt"
Yfilename = "LastY.txt"
GifOutputName = "Dynamics.gif"

# Is warm-start data available?
WarmStart = true

# Add noise to prey data?
ShakeUpX = false

# Set a grid size for simulation
GridSize = 512

# Choose the time length of the simulation
TMAX = 24

# Choose the time point to begin animating from (i.e., animate TSTART:TMAX; note min is 1, not 0)
TSTART = 1

################################################################################
### PARAMETERS ###

μ = 4
β = 5
D1 = 0.001
D2 = 0.2

δ1x = 0 # prey min initialization
δ2x = 1 # prey max initialization
δ1y = 0.3 # predator min initialization
δ2y = 0.4 # predator max initialization

################################################################################
### SET UP GRID ###

# Actual grid size
ShadowGridSize = GridSize+2 # +2 because julia is 1-indexed; the size of the grid including the PBC
GridEnd = GridSize+1 # +1 because julia is 1-indexed; the grid end point excluding the PBC

# Shadow grid with PBC elements
I1S = collect(1:ShadowGridSize) # Ghost elements for grid PBC
I2S = I1S # alias for second dimension

# True grid without PBC elements
I1 = collect(2:GridEnd) # This is one of those times 1-indexing is annoying; the grid points excluding the PBC
I2 = I1 # alias for second dimension

# Number of time points
T = collect(1:TMAX)
TANIM = collect(TSTART:TMAX)

# time points where a progress message is displayed
TPrint = []
TPrintAnim = []
UpdateInterval = 10
for i = 1:UpdateInterval
    push!(TPrint,Int(floor((i/UpdateInterval)*TMAX)))
    push!(TPrintAnim,Int(floor((i/UpdateInterval)*(TMAX-TSTART)+TSTART)))
end

################################################################################
### VARIABLES ###

X = zeros(ShadowGridSize,ShadowGridSize,TMAX)
Y = zeros(ShadowGridSize,ShadowGridSize,TMAX)
∇X = zeros(ShadowGridSize,ShadowGridSize,TMAX)
∇Y = zeros(ShadowGridSize,ShadowGridSize,TMAX)

################################################################################
### SIMULATION FUNCTIONS ###

function x_update(t)
    """
    Prey dynamics
    """
    for i1 in I1
        for i2 in I2
            X[i1,i2,t] = μ*X[i1,i2,t-1]*(1 - X[i1,i2,t-1])*exp(-β*Y[i1,i2,t-1]) + D1*∇X[i1,i2,t-1]
            if X[i1,i2,t] < 0
                X[i1,i2,t] = 0
            end
        end
    end
end

function y_update(t)
    """
    Predator dynamics
    """
    for i1 in I1
        for i2 in I2
            Y[i1,i2,t] = X[i1,i2,t-1]*(1 - exp(-β*Y[i1,i2,t-1])) + D2*∇Y[i1,i2,t-1]
            if Y[i1,i2,t] < 0
                Y[i1,i2,t] = 0
            end
        end
    end
end

function delx_update(t)
    """
    Spatial coupling
    """
    for i1 in I1
        for i2 in I2
            ∇X[i1,i2,t] = X[i1-1,i2,t] + X[i1+1,i2,t] + X[i1,i2-1,t] + X[i1,i2+1,t] - 4*X[i1,i2,t]
        end
    end
end

function dely_update(t)
    """
    Spatial coupling
    """
    for i1 in I1
        for i2 in I2
            ∇Y[i1,i2,t] = Y[i1-1,i2,t] + Y[i1+1,i2,t] + Y[i1,i2-1,t] + Y[i1,i2+1,t] - 4*Y[i1,i2,t]
        end
    end
end

function PBC(t)
    """
    Period Boundary Condition
    """
    for i1 in I1
        X[i1,1,t] = X[i1,GridEnd,t]
        X[i1,ShadowGridSize,t] = X[i1,2,t]
        Y[i1,1,t] = Y[i1,GridEnd,t]
        Y[i1,ShadowGridSize,t] = Y[i1,2,t]
    end
    for i2 in I2
        X[1,i2,t] = X[GridEnd,i2,t]
        X[ShadowGridSize,i2,t] = X[2,i2,t]
        Y[1,i2,t] = Y[GridEnd,i2,t]
        Y[ShadowGridSize,i2,t] = Y[2,i2,t]
    end
end

function XAddNoise(t)
    """
    Add noise to data; i.e., if importing results that are steady-state.
    Noise centred at 0, uniform in +/- δ2x (upper initialization limit)
    After noise addition, values >1 -> 1, value <0 -> 0
    """
    X[:,:,t] += 2*(rand(ShadowGridSize,ShadowGridSize)-0.5*ones(ShadowGridSize,ShadowGridSize))*δ2x
    for i1 in I1
        for i2 in I2
            if X[i1,i2,t] > 1
                X[i1,i2,t] = 1
            elseif X[i1,i2,t] < 0
                X[i1,i2,t] = 0
            end
        end
    end
end

################################################################################
### INITIALIZE FOR t=1 ###
if WarmStart
    # Assign warm start data as initial value
    X[:,:,1] = readdlm(Xfilename, '|', Any)
    Y[:,:,1] = readdlm(Yfilename, '|', Any)
else
    # Assign random values
    X[:,:,1] = rand(ShadowGridSize,ShadowGridSize)*(δ2x-δ1x) + ones(ShadowGridSize,ShadowGridSize)*δ1x
    Y[:,:,1] = rand(ShadowGridSize,ShadowGridSize)*(δ2y-δ1y) + ones(ShadowGridSize,ShadowGridSize)*δ1y
end

# Add noise to prey data; useful to add disturbance to steady-state warm-start results
if ShakeUpX
    XAddNoise(1)
end

################################################################################
### RUN SIMULATION ###
println("Starting simulation...")

# Iterate over time horizon starting from t=1
for t = 2:TMAX
    # Apply boundary conditions to previous time step
    PBC(t-1)

    # Calculate derivative analogues (using PBC)
    delx_update(t-1)
    dely_update(t-1)

    # Update X and Y for period t
    x_update(t)
    y_update(t)

    if t in TPrint
        println("Simulation progress: "*string(Int(floor(100*t/TMAX)))*"%")
    end
end
println("Done simulation!")

# Save last time point for warm start
file = open(Xfilename,"w")
    writedlm(file, X[:,:,TMAX], '|')
close(file)
file = open(Yfilename,"w")
    writedlm(file, Y[:,:,TMAX], '|')
close(file)

################################################################################
### ANIMATE ###
# Start prey animation loop
println("Animating dynamics...")
anim = @animate for t in TANIM
    # Takes the t-th slice and plot
    # for a different colorscheme: color=:halide or, color=:magma
    h1 = heatmap(X[:,:,t], clims=(0,1), xlabel="Prey")
    h2 = heatmap(Y[:,:,t], clims=(0,1), xlabel="Predators")
    plot(h1,h2, layout=(1,2), legend=false, wsize=(1200,600))

    # Add the number to the heatmap
    annotate!(0,0,string(lpad(t - (TSTART-1),ndigits(TMAX),"0")))

    # Progress message:
    if t in TPrintAnim
        println("Animation progress: "*string(Int(floor(100*(t-TSTART)/(TMAX-(TSTART-1)))))*"%")
    end
end
gif(anim, GifOutputName, fps = 7)
println("Done animation!")
