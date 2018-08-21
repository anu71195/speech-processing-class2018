window=10000;
threshold=4000;
iteration_no=0;
file=open("0_1_10iterations_without_charging.txt","r")
filename="150101010";
print (file)
indexx=0;
start=0;
count=0;
end=start+window;
complete=0;
data=[];
flag=0;
for line in file:
	data.append(int(line.split('\n')[0]));
while(1):
	for i in range(start,len(data)+1):
		if(data[i]>threshold):
			first_high=i;
			last_high=i;
			start=first_high+1;
			break;
	while (1):
		if(start+window >=len(data)):
			complete=1;
			break;
		for i in range(start,start+window):
			if(data[i]>threshold):
				last_high=i;
				flag=1;
		start=last_high+1;
		if(flag==0):
			break;
		flag=0;
	if(count%10==0):
		iteration_no+=1;
	count+=1;
	if(count>=101):
		break;
	fw=open(filename+"_"+str((count-1)%10)+"_"+str(iteration_no)+".txt","w")
	for i in range(first_high-int(window/5),last_high+int(window/5)):
		fw.write(str(data[i]));
		fw.write("\n");
	if(complete):
		break;

	

