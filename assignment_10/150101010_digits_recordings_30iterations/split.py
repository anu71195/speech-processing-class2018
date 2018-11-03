window=10000;
threshold=2000;
iteration_no=0;
file=open("digits.txt","r")
filename="150101010";
print (file)
indexx=0;
start=0;
count=0;
my_count=1;
end=start+window;
complete=0;
data=[];
flag=0;
digits=["0","1","2","3","4","5","6","7","8","9"];

for line in file:
	data.append(int(line.split('\n')[0]));
print(len(data))
while(1):
	print(start,my_count);
	my_count+=1;
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
	#if(count>=51):
#		break;
	fw=open(filename+"_"+digits[(count-1)%10]+"_"+str(iteration_no)+".txt","w")
	print(filename+"_"+digits[(count-1)%10]+"_"+str(iteration_no)+".txt","w")
	for i in range(first_high-int(window/5),last_high+int(window/5)):
		fw.write(str(data[i]));
		fw.write("\n");
	if(complete):
		break;

	

