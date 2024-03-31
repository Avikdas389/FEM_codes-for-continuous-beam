# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 12:11:17 2024

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

user_input = input("Enter length of elements in meter separated by spaces: ")
l = user_input.split()
l = [float(element) for element in l]
print("Length of elements are:", l)
N_element= len(l)
N_node = N_element+1
E = float(input("Enter the Elastic Modulus in Pa :"))
#E= 2e11
I = float(input("Enter the Moment of Inertia m^4 :"))
#I = 0.000004
EI = E*I
N_div=17
T=np.zeros(N_node)
for j in range(0,N_element):
    T[j+1] = l[j] + T[j]
K=np.zeros([2*N_node,2*N_node])
for i in range(N_element):
    n = i
    f = n+1
    K[2*n,2*n]=K[2*n,2*n]+12/(l[i]**3)
    K[2*n,2*n+1]=K[2*n,2*n+1]+6/(l[i]**2)
    K[2*n,2*f]=K[2*n,2*f]-12/(l[i]**3)
    K[2*n,2*f+1]=K[2*n,2*f+1]+6/(l[i]**2)
    
    K[2*n+1,2*n]=K[2*n+1,2*n]+6/(l[i]**2)
    K[2*n+1,2*n+1]=K[2*n+1,2*n+1]+4/l[i]
    K[2*n+1,2*f]=K[2*n+1,2*f]-6/(l[i]**2)
    K[2*n+1,2*f+1]=K[2*n+1,2*f+1]+2/l[i]
    
    K[2*f,2*n]=K[2*f,2*n]-12/(l[i]**3)
    K[2*f,2*n+1]=K[2*f,2*n+1]-6/(l[i]**2)
    K[2*f,2*f]=K[2*f,2*f]+12/(l[i]**3)
    K[2*f,2*f+1]=K[2*f,2*f+1]-6/(l[i]**2)
    
    K[2*f+1,2*n]=K[2*f+1,2*n]+6/(l[i]**2)
    K[2*f+1,2*n+1]=K[2*f+1,2*n+1]+2/l[i]
    K[2*f+1,2*f]=K[2*f+1,2*f]-6/(l[i]**2)
    K[2*f+1,2*f+1]=K[2*f+1,2*f+1]+4/l[i]
          
print('Global stiffness matrix')    
print(K*E*I)



#******Loading ******#
Tr=np.zeros(N_element)
P= np.zeros(N_div*N_element)
y_2= np.zeros(N_div)
F=np.zeros(2*N_node)
V_end=np.zeros(2*N_node-2)
M_fend=np.zeros(2*N_node-2)


for j in range(N_element):
   
       t = int(input(f"Loading pattern for element {j+1} type 1 only for UDL:"))             
       if (t == 1):
           u= int(input("Enter the UDL value in N/m :"))
           Tr[j]= 1 
           
           F[2*j]= F[2*j] - u*l[j]/2
           F[2*j+1]= F[2*j+1] - u*l[j]**2/12
           F[2*j+2]= F[2*j+2] - u*l[j]/2
           F[2*j+3]= F[2*j+3] + u*l[j]**2/12
           
           M_fend[2*j]= M_fend[2*j] - u*l[j]**2/12 
           M_fend[2*j+1]= M_fend[2*j+1] + u*l[j]**2/12
           
           V_end[2*j]= V_end[2*j]+ u*l[j]/2
           V_end[2*j+1]= V_end[2*j+1]- u*l[j]/2
           
           x_2= np.linspace(0,l[j],N_div)
           for i in range(0,N_div):
               y_2[i]= u*l[j]*x_2[i]/2  - u*x_2[i]**2/2
           g=0
           for q in range(N_div*j,N_div*j+N_div):
               P[q]= y_2[g]
               g= g +1  
       else :
           p= int(input("Enter the Concentrated load value acting at midpoint in N:"))
          
           F[2*j]= F[2*j] - p/2
           F[2*j+1]= F[2*j+1] - p*l[j]/8
           F[2*j+2]= F[2*j+2] - p/2
           F[2*j+3]= F[2*j+3] + p*l[j]/8
           
           M_fend[2*j]= M_fend[2*j] - p*l[j]/8
           M_fend[2*j+1]= M_fend[2*j+1] + p*l[j]/8
           
           V_end[2*j]= V_end[2*j] + p/2
           V_end[2*j+1]= V_end[2*j+1] - p/2
          
           x_2= np.linspace(0,l[j],N_div)
           for i in range(0,N_div):
               if i < N_div/2:
                   y_2[i]= p*x_2[i]/2
               else:
                   y_2[i]= p*l[j]/2 - p*x_2[i]/2
           
           g=0
           for q in range(N_div*j,N_div*j+N_div):
               P[q]= y_2[g]
               g= g +1

print('Global load vector')
print(F)


#****Giving settlement as input****/

user_se = input("Enter the node numbers where settlement is provided: ")
se = user_se.split()
se = [int(element) for element in se]
N_set = len(se)

K= K*E*I
user_inp = input("Enter the value of settlement sequentially in m : ")
Set = user_inp.split()
Set = [float(element) for element in Set]
for i in range(N_set):
    o = se[i]
    K[2*o,2*o]= K[2*o,2*o] + 10e9
    F[2*o]= F[2*o]+ 10e9*Set[i]
    
print("Updated Stiffnex matrix with Penalty Approach" )   
print(K)


#******Boundary Conditions ********#


user_inpt = input("Enter Degree of freedom where Displacement is zero sequentially : ")
d = user_inpt.split()
d = [int(element) for element in d]
N_d= len(d)

for i in range(N_d):
    s= d[i]-i
    
    K=np.delete(K, s, 0)  
    K=np.delete(K, s, 1)
    F=np.delete(F,s,0)



print('Stiffness Matrix after elimination')
print (K)

#******Solutions******#

print('Global load vector after elimination')
print(F)
All_DOF= np.zeros(2*N_node)
Q= np.linalg.solve(K, F)
Qg=np.linspace(0,2*N_node-1,2*N_node) # all the dof indices
All_DOF=np.linspace(0,0,2*N_node)# displacements at all the dofs
NZQ=np.setdiff1d(Qg,d)  # the array containing the indices of the free dofs
ZQ=np.setdiff1d(Qg,NZQ)  # the array containing the indices of the constrained dofs

nNZDOF= 2*N_node - N_d

for i in range(nNZDOF):
  j=int(NZQ[i])
  All_DOF[j]=Q[i]
print('Displacement at all dofs')
print(All_DOF)

#**** Calculation of Moment *******/
All_DOF=-All_DOF # As sign convention is different for moment Calculation
M_final=np.zeros(2*N_node-2)
for i in range(0,N_element):
  for j in range(i*2,i*2+2):
     if j % 2 == 0:
      M_final[j]= M_fend[j] + 4*E*I*All_DOF[j+1]/l[i] + 2*E*I*All_DOF[j+3]/l[i] -6*E*I*(All_DOF[j+2]-All_DOF[j])/l[i]**2
     else:
      M_final[j]= M_fend[j] + 4*E*I*All_DOF[j+2]/l[i] + 2*E*I*All_DOF[j]/l[i] -6*E*I*(All_DOF[j+1]-All_DOF[j-1])/l[i]**2  

print("The Fixed end Moments are")
print(M_fend)
print(" The Final end Moments are")
print(M_final)

#******* Plotting of BMD individuals*****/
M_plot= np.zeros(N_node)
M_plot[0]= M_final[0]
M_plot[N_node-1]= -M_final[2*N_node-3]
z=1
u1=0


for i in range(1,2*N_node-3):   
    if i % 2 == 0:
        u1=u1+1
    else:
        M_plot[z]= -M_final[i]
        z=z+1

print("Moments considered for plotting ")
print(M_plot)
P_z= np.zeros(N_div)

for i in range(0,N_element):
    x_values= np.linspace(T[i],T[i+1], N_div)
    y_1= np.linspace(M_plot[i],M_plot[i+1],N_div)
    plt.plot(x_values,y_1, linestyle='-', color = 'blue',label= ' Variation of End Moments ')
    plt.xlabel('Variation along length (m)')
    plt.ylabel('Moment (N-m)')
    plt.title('Bending Moment Diagram')
    plt.grid(True)  # Add grid

for i in range(0,N_element):
    x_values= np.linspace(T[i],T[i+1], N_div)
    P_z = np.zeros(N_div)  # Initialize P_z array for each iteration
    for j in range(N_div):
        P_z[j] = P[N_div*i + j]        
    plt.plot(x_values,P_z, marker='o',linestyle='-', color = 'red')

#***** Plotting of the central line //
z_1= np.zeros(len(T))
plt.plot(T,z_1,linestyle='-', color = 'green')
plt.legend(['Variation of End Moments'])
plt.show()

#******* Plotting of BMD Total***/  
for i in range(0,N_element):
    x_values= np.linspace(T[i],T[i+1], N_div)
    for j in range(N_div):
        P_z[j] = P[N_div*i + j]
    y_1= np.linspace(M_plot[i],M_plot[i+1],N_div)
    if i % 2 == 0:
       plt.plot(x_values,P_z +y_1, linestyle='-', color = 'black')
    else:
       plt.plot(x_values,P_z +y_1, linestyle='-', color = 'black')    
    plt.xlabel('Variation along length (m)')
    plt.ylabel('Moment (N-m)')
    plt.title('Bending Moment Diagram ')
    plt.grid(True)

#***** Plotting of the central line //
z_1= np.zeros(len(T))
plt.plot(T,z_1,linestyle='-', color = 'green')

plt.show()

#***** Plotting of SFD *** //
V_pl1= np.zeros(4)
M_s= np.zeros(N_element)

for i in range(0,N_element):
  for j in range(i*2,i*2+1):
      M_s[i]= M_final[j]/l[i] + M_final[j+1]/l[i]

for i in range(0,N_element):
    for j in range(i*2,i*2+2):
      V_end[j]= V_end[j] - M_s[i]

print("Shear Force at the joints " )
print(V_end)
#print(Tr)
#joining of end points #
p01 = [T[0], T[0]]
q01 = [V_end[0], 0] # y-cordinate
plt.plot(p01, q01, 'red', linestyle="-")
p01 = [T[N_node-1], T[N_node-1]]
q01 = [V_end[2*N_element-1], 0]
plt.plot(p01, q01, 'red', linestyle="-")


for j in range (1,N_element):
    p01 = [T[j], T[j]]
    q01 = [V_end[2*j-1], V_end[2*j]]
    plt.plot(p01, q01, 'red', linestyle="-") 
       
for j in range(N_element):
    if Tr[j] == 0 :
        V_pl1[0]= V_end[2*j]
        V_pl1[1]= V_end[2*j]
        V_pl1[2]= V_end[2*j+1]
        V_pl1[3]= V_end[2*j+1]
        
        x_1= np.linspace(T[j],T[j+1],3)
        p01 = [x_1[0], x_1[1]]
        q01 = [V_pl1[0], V_pl1[1]]
        plt.plot(p01, q01, 'red', linestyle="-")
        p01 = [x_1[1], x_1[1]]
        q01 = [V_pl1[1], V_pl1[2]]
        plt.plot(p01, q01, 'red', linestyle="-")
        p01 = [x_1[1], x_1[2]]
        q01 = [V_pl1[2], V_pl1[3]]
        plt.plot(p01, q01, 'red', linestyle="-")
        
    else :
        
        p01 = [T[j], T[j+1]]
        q01 = [V_end[2*j], V_end[2*j+1]]
        plt.plot(p01, q01, 'red', linestyle="-")    
    plt.xlabel('Variation along length (m)')
    plt.ylabel('Shear Force (N)')
    plt.title('Shear Force Diagram')
    plt.grid(True)

z_1= np.zeros(len(T))
plt.plot(T,z_1,linestyle='-', color = 'green')        

plt.show()

