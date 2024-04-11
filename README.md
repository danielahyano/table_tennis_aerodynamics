# Table Tennis Aerodynamics Project 

The trajectory of table tennis balls is affected by the gravitational force, lift and drag. In professional table tennis, backspin and topspin are oftentimes combined with sidespin, producing interesting effects on the trajectory of the ball. One of the most important pieces of the game are the serves, where the player who is serving has the freedom to inject the desired amount of spin and velocity to the ball. In this project, I plan to study the aerodynamics of table tennis serves, by simulating the trajectory of the table tennis balls, taking into account the different forces and comparing the trajectory for different initial conditions.

Source: https://pubs.aip.org/aapt/ajp/article/88/11/934/1058353/Flight-and-bounce-of-spinning-sports-balls

Milestones:
- [X] Simulation of flight in 2D (without spin)
- [X] Simulation of flight in 2D (with topspin and backspin)
- [X] Simulation of flight trajectory in 3D
- [X] Simulation of flight trajectory in 3D with sidespin, topspin and backspin
- [ ] Make modifications to the bounce, taking into account the dissipation of energy in form of heat and sound. 

Steps to run the program:
1. Compile:
'g++ -o myprog balle.cpp'
2. Run the program:
'./myprog '
3. To plot you can use two python programs:
'plot_tt_animated.py' which is an animation 
and
'plot_tt.py' which is a simple plot.
