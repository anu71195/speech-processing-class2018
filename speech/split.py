file=open("0_1_10iterations_without_charging.txt","r")
filename="iteration";
print (file)
indexx=0;
start=0;
count=0;
window=27000;
end=start+window;
threshold=4000;
data=[];
flag=0;
for line in file:
	data.append(int(line.split('\n')[0]));
while(1)	:
	start=end;
	end=start+window;
	if(end>=len(data)-1):
		break;
	for i in range(start,end+1):
		if(data[i]>threshold):
			flag=1;
			count+=1;
			fw=open(filename+str(count),"w")
			for j in range(start,end+1):
				fw.write(str(data[j]))
				fw.write("\n");
			break;

	# if (flag):
		# flag=0;
		# continue;
	# count+=1;
	
print(count)