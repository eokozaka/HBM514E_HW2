double sort(double *a, int n){
	bool sorted=false,flag = false;
	int temp=0,swap=0;
	while(sorted==false){
		flag = false;
		for(int i=0;i<n;i++){
			if(a[i]<a[i+1]){
				flag = true;
				temp = a[i];
				a[i] = a[i+1];
				a[i+1] = temp;
			}
		}
		if(flag == false){sorted = true;}
	}
return *a;
}
