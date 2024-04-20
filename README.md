# Table Tennis Aerodynamics Project 

The trajectory of table tennis balls is affected by the gravitational force, lift and drag. In professional table tennis, backspin and topspin are oftentimes combined with sidespin, producing interesting effects on the trajectory of the ball. One of the most important pieces of the game are the serves, where the player who is serving has the freedom to inject the desired amount of spin and velocity to the ball. In this project, I plan to study the aerodynamics of table tennis serves, by simulating the trajectory of the table tennis balls, taking into account the different forces and comparing the trajectory for different initial conditions.

I have based my project on the paper: https://pubs.aip.org/aapt/ajp/article/88/11/934/1058353/Flight-and-bounce-of-spinning-sports-balls

I have used as the starting point of ny balle.cpp code from http://www.algarcia.org/nummeth/nummeth.html. 

Milestones:
- [X] Simulation of flight in 2D (without spin)
- [X] Simulation of flight in 2D (with topspin and backspin)
- [X] Simulation of flight trajectory in 3D
- [X] Simulation of flight trajectory in 3D with sidespin, topspin and backspin
- [ ] Make modifications to the bounce, taking into account the dissipation of energy in form of heat and sound. 

- Next steps:
 - Adapt so that you can run the code for different sports (just changing parameters, but making it automated)
 - Compare results with data to check if code looks good

Steps to run the program:
1. Compile:
'g++ -o myprog balle.cpp'
2. Run the program:
'./myprog '
3. To plot you can use one of my three python programs:
'plot_tt_animated.py' which is an animation that shows the plot in z-x axis
'plot_tt_up_vision.py' which is an animation from the "up" vision
and 'plot_tt.py' which is a simple plot to see the trajectory!.
