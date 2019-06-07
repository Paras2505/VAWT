MIN_Line=3

cfd_file= open('forces.dat','r')
list=cfd_file.readlines()
# k=list[3].split(')')
# l=k[4].split()[2]
# print(l)
MAX_Line=len(list)
cost=0
if(len(list)>=MAX_Line):
    for i in range(MIN_Line,MAX_Line):
        print(float(list[i].split(')')[4].split()[2]),' ',float(list[i].split(')')[5].split()[2]))
        cost+=float(list[i].split(')')[4].split()[2])+float(list[i].split(')')[5].split()[2])
    cost=cost/(MAX_Line-MIN_Line)

print(cost)
cfd_file.close()
