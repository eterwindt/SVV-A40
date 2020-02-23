import numpy as np
from numpy import cos,sin,tan,pi,sqrt, arcsin



#programm starts here
Ca=505
ha=161
hst=13
tsk=1.1
tsp=2.4
tst=1.2
wst=17
nst=11
listz=[]
listy=[]
#Circumference of aileron cross-section
O_ail=pi*0.5*ha+2*sqrt((Ca-0.5*ha)**2+(0.5*ha)**2)
#calculating the stringer pitch
str_pitch=(O_ail)/(nst)

#Calculating coordinates of stringers on semi-circle (only for positive y)
theta_str_pitch=2*str_pitch/ha
    #the coordiates of the stringers 
for i in range(nst):
    if theta_str_pitch*i<(0.5*pi):
        z_str=0.5*ha*cos(theta_str_pitch*i)-0.5*ha
        y_str=0.5*ha*sin(theta_str_pitch*i)
        listz.append(z_str)
        listy.append(y_str)
#the remaining length in the semi circle
length_semi_remain=0.25*ha*pi-(len(listz)-1)*str_pitch

#Calculation for the coordinates of the stringers on slope side
gamma=np.arctan((0.5*ha)/(Ca-0.5*ha))
#next stringer on slope after semi circle
z_str=-1.0*cos(gamma)*(str_pitch-length_semi_remain)-0.5*ha
y_str=0.5*ha-sin(gamma)*(str_pitch-length_semi_remain)
listz.append(z_str)
listy.append(y_str)
#area stringer
A_str=(wst+hst)*tst


#remaining stringers on slope
for i in range((int(nst/2)+1)-len(listz)):
    z_str=z_str-cos(gamma)*str_pitch
    listz.append(z_str)
for i in range((int(nst/2)+1)-len(listy)):
    y_str=y_str-sin(gamma)*str_pitch
    listy.append(y_str)

listz.extend(np.flip(listz[1:7]))
listy.insert(2,0.5*ha)
listy.extend(-np.flip(listy[1:7]))



My=500
Mz=200
Iyy=1000
Izz=2000
siglist=[]
#for i in range(12):
#    sig=listz[i]*My/Iyy-listy[i]*Mz/Izz
#    siglist.append(sig)

Az_str=[A_str*i for i in listz]

A_sstr=11*A_str
A_sp=tsp*ha
A_wbc=0.5*((0.5*ha)**2*pi-(0.5*ha-tsk)**2*pi)
A_wbs=2*tsk*((Ca-0.5*ha)**2+(0.5*ha)**2)**0.5
As=A_sstr+A_sp+A_wbc+A_wbs

Az_sp=-A_sp*0.5*ha
Az_wbc=A_wbc*(-2*0.5*ha/pi)
Az_wbs=-A_wbs*(0.5*ha+0.5*Ca)
Azs=sum(Az_str)+Az_sp+Az_wbc+Az_wbs            
z=Azs/As


listz.insert(2,-0.5*ha)
listz.insert(11,-0.5*ha)

# Coordinate to theta conversion
angle_list = []
for n in range(0,3):
    z_coor = -listz[n]-0.5*ha
    if (z_coor) < 0:
        angle = np.arctan((listy[n]) / (z_coor) ) + pi
    else:
        angle = 0.5*pi
    angle_list.append(angle)
angle_list.insert(0,1.5*pi)
angle_list.insert(1,2*pi-angle_list[2])

#boom areas: Creating a list of boom areas for the cases My and Mz seperately
list_str_pitch=[str_pitch]*7
list_str_pitch.insert(0,str_pitch-length_semi_remain)
list_str_pitch.insert(8,str_pitch-length_semi_remain)
listz=[i-z for i in listz]
b=list_str_pitch
B_list_Mz=[]
B_list_My=[]
for i in range(8):
    B_Mz=A_str+tsk*b[i]/6*(2+listy[i+2]/listy[i+3])+tsk*b[i+1]/6*(2+listy[i+4]/listy[i+3])
    B_list_Mz.append(B_Mz)
    B_My=A_str+tsk*b[i]/6*(2+listz[i+2]/listz[i+3])+tsk*b[i+1]/6*(2+listz[i+4]/listz[i+3])
    B_list_My.append(B_My)
B_Mz_sp=tsk*b[0]/6*(2+listy[3]/listy[2])+tsp*ha/6*(2+listy[11]/listy[2]) #spar
B_My_sp=tsk*b[0]/6*(2+listz[3]/listz[2])+tsp*ha/6*(2+listz[11]/listz[2]) #spar
B_list_Mz.insert(0,B_Mz_sp) #adding sparbooms to the lists
B_list_Mz.insert(9,B_Mz_sp)
B_list_My.insert(0,B_My_sp)
B_list_My.insert(9,B_My_sp)

B1_Mz=[]
B2_Mz=[]
B1_My=[]
B2_My=[]
#for i in range(4):
#    B1_area_Mz=(0.5*ha*tsk*(0.5*(theta2-theta1)+0.25*(np.sin(theta1*2)-np.sin(theta2*2))+(np.cos(theta2)-np.cos(theta1))*(np.sin(theta2))))/((np.sin(theta1))**2-np.sin(theta1)*np.sin(theta2))
#    B2_area_Mz=(0.5*ha*tsk*(0.5*(theta2-theta1)+0.25*(np.sin(theta1*2)-np.sin(theta2*2))+(np.cos(theta2)-np.cos(theta1))*(np.sin(theta1))))/((np.sin(theta2))**2-np.sin(theta1)*np.sin(theta2))

#B1_area_My=(0.5*ha*tsk*((0.5*ha)**2*(0.5*(theta2-theta1)+0.25*(sin(2*theta2)-sin(2*theta1)))-2*X_c*0.5*ha*(sin(theta2)-sin(theta1))+(X_c)**2*(theta2-theta1))-0.5*ha*tsk*(0.5*ha*(sin(theta2)-sin(theta1))+X_c*(theta1-theta2))*(cos(theta2)*0.5*ha-X_c))/((cos(theta1)*0.5*ha-X_c)**2-(cos(theta1)*0.5*ha-X_c)*(cos(theta2)*0.5*ha-X_c))
#B2_area_My=(0.5*ha*tsk*((0.5*ha)**2*(0.5*(theta2-theta1)+0.25*(sin(2*theta2)-sin(2*theta1)))-2*X_c*0.5*ha*(sin(theta2)-sin(theta1))+(X_c)**2*(theta2-theta1))-0.5*ha*tsk*(0.5*ha*(sin(theta2)-sin(theta1))+X_c*(theta1-theta2))*(cos(theta1)*0.5*ha-X_c))/((cos(theta2)*0.5*ha-X_c)**2-(cos(theta1)*0.5*ha-X_c)*(cos(theta2)*0.5*ha-X_c))

#Add function for boom areas of circular section here!!!!
B_list_Mz.insert(0,1) #temporary value for circular section
B_list_Mz.insert(0,1) #temporary calue for circular section
B_list_Mz.insert(13,1)#temporary calue for circular section
B_list_My.insert(0,1) #temporary calue for circular section
B_list_My.insert(0,1) #temporary calue for circular section
B_list_My.insert(13,1)#temporary calue for circular section

#moment of inertias
listy2=[i ** 2 for i in listy]
listz2=[i ** 2 for i in listz]
Izz=np.dot(listy2,B_list_Mz)
Iyy=np.dot(listz2,B_list_My)
V_y=-15000
V_z=0
#q_b_Mz=[-V_y/Izz*B_list_Mz[i]*listy)

#contribution van de semi-circle naar de boom area's 
#eerst angles (theta) van alle punten om de semi-circle berekenen
#Met afwijkende nummering van booms (boom12,boom13,boom1,boom2,boom3)

angle_lijst=[1.5*pi,4.40,pi,1.88,0.5*pi]
#Skin contribution of semi-circle of Mz idealization

# for Mz (Boom12,Boom13,Boom1,Boom2,Boom3)
semiskcontri_Mz=[]

for i in range(5):
    theta2=angle_lijst[i]
    if i<3:
        theta1=angle_lijst[i+1]
        B2_area_Mz=(0.5*ha*tsk*(0.5*(theta2-theta1)+0.25*(np.sin(theta1*2)-np.sin(theta2*2))+(np.cos(theta2)-np.cos(theta1))*(np.sin(theta1))))/((np.sin(theta2))**2-np.sin(theta1)*np.sin(theta2))
    else:
        B2_area_Mz=0.0
    theta1=angle_lijst[i]
    if i>0:
        theta2=angle_lijst[i-1]
        B1_area_Mz=(0.5*ha*tsk*(0.5*(theta2-theta1)+0.25*(np.sin(theta1*2)-np.sin(theta2*2))+(np.cos(theta2)-np.cos(theta1))*(np.sin(theta2))))/((np.sin(theta1))**2-np.sin(theta1)*np.sin(theta2))
    else:
        B1_area_Mz=0.0
    
    jup=B2_area_Mz+B1_area_Mz
    semiskcontri_Mz.append(jup)
semiskcontri_Mz[2]=0.0
# for My (Boom12,Boom13,Boom1,Boom2,Boom3)
semiskcontri_My=[]
X_c=listz[0]-0.5*ha
for i in range(5):
    theta2=angle_lijst[i]
    if i<3:
        theta1=angle_lijst[i+1]
        B2_area_My=(0.5*ha*tsk*((0.5*ha)**2*(0.5*(theta2-theta1)+0.25*(sin(2*theta2)-sin(2*theta1)))-2*X_c*0.5*ha*(sin(theta2)-sin(theta1))+(X_c)**2*(theta2-theta1))-0.5*ha*tsk*(0.5*ha*(sin(theta2)-sin(theta1))+X_c*(theta1-theta2))*(cos(theta1)*0.5*ha-X_c))/((cos(theta2)*0.5*ha-X_c)**2-(cos(theta1)*0.5*ha-X_c)*(cos(theta2)*0.5*ha-X_c))
    else:
        B2_area_My=0.0
    theta1=angle_lijst[i]
    if i>0:
        theta2=angle_lijst[i-1]
        B1_area_My=(0.5*ha*tsk*((0.5*ha)**2*(0.5*(theta2-theta1)+0.25*(sin(2*theta2)-sin(2*theta1)))-2*X_c*0.5*ha*(sin(theta2)-sin(theta1))+(X_c)**2*(theta2-theta1))-0.5*ha*tsk*(0.5*ha*(sin(theta2)-sin(theta1))+X_c*(theta1-theta2))*(cos(theta2)*0.5*ha-X_c))/((cos(theta1)*0.5*ha-X_c)**2-(cos(theta1)*0.5*ha-X_c)*(cos(theta2)*0.5*ha-X_c))
    else:
        B1_area_My=0.0
    
    jup=B2_area_My+B1_area_My
    semiskcontri_My.append(jup)
    
    
    
#adding the semi circle skin contribution to the boom area lists
B_list_Mz[0]=semiskcontri_Mz[2]
B_list_Mz[1]=semiskcontri_Mz[3]
B_list_Mz[2]=semiskcontri_Mz[4]+B_list_Mz[2]
B_list_Mz[11]=semiskcontri_Mz[0]+B_list_Mz[11]
B_list_Mz[12]=semiskcontri_Mz[1]

B_list_My[0]=semiskcontri_My[2]
B_list_My[1]=semiskcontri_My[3]
B_list_My[2]=semiskcontri_My[4]+B_list_My[2]
B_list_My[11]=semiskcontri_My[0]+B_list_My[11]
B_list_My[12]=semiskcontri_My[1]

#adding stiffener area to boom area's
A_str=(wst+hst)*tst

for i in range(12):
    B_list_Mz[i]=B_list_Mz[i]+A_str
    B_list_My[i]=B_list_My[i]+A_str

B_list_Mz[0]=0.0 #Ask TA about this

#Calculation of Izz and Iyy
list_Izz=[]
list_Iyy=[]
for i in range(12):
    Izz=B_list_Mz[i]*(listy[i])**2
    list_Izz.append(Izz)
IZZ=sum(list_Izz)

for i in range(12):
    Iyy=B_list_My[i]*(listz[i])**2
    list_Iyy.append(Iyy)
IYY=sum(list_Iyy)