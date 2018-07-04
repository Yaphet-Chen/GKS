The attached is a Very First Version of a light 1d high order GKS program, contains:

WENO5Z,
KFVS1st, KFVS2nd, GKS1st, GKS2nd, HLLC,
RK4 and multistage method(Stage2Oorder4 ),
OMP parallel






g++ -o3 -std=c++0x -fopenmp -o gks1d highordergks1d.cpp