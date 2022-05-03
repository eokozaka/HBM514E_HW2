double sort(int *a, int n){
	bool sorted=false;
	int temp=0,swap=0;
	while(sorted==false){
		for(int i=0;i<n;i++){
			if(a[i]<a[i+1]){
				temp = a[i];
				a[i] = a[i+1];
				a[i+1] = temp;
			}

		}
	}
return *a;
}
