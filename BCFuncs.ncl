undef("ran1")
function ran1(idum[1]:integer)
local IA,IM,AM,IQ,IR,NTAB,NDIV,EPS,RNMX,iv,iy,j,k,output
begin
	IA = 16807
	IM = 2147483647
	AM = 1./IM
	IQ = 127773
	IR = 2836
	NTAB = 32
	NDIV = 1+(IM-1)/NTAB
	EPS = 1.2e-7
	RNMX = 1.-EPS
	iy = 0
	iv = new(NTAB,integer)
	
	if ((idum.le.0).or.(.not.iy)) then
		if (-(idum).lt.1) then
			idum = 1
		else
			idum = -(idum)
		end if
		do j = NTAB+7,0,1
			k = floattoint(idum/IQ)
			idum = IA*(idum-k*IQ)-IR*k
			if (idum.lt.0) then idum = idum+IM end if
			if (j.lt.NTAB) then iv(j) = idum end if
		end do
		iy = iv(0)
	end if
	k = idum/IQ
	idum = IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) then idum = idum + IM end if
	j = iy/NDIV
	iy = iv(j)
	iv(j) = idum
	if (AM*iy.gt.RNMX)
		output = RNMX
	else
		output = AM*iy
	end if
	return(output)
end
	
undef("gasdev")
function gasdev(idum[1]:integer)
local v1,v2,r,gset,iset,fac,output
begin
	iset = 0
	if (iset.eq.0) then
		v1 = 2.*ran1(idum)-1.
		v2 = 2.*ran1(idum)-1.
		r = v1*v1 + v2*v2
		do while (r.ge.1.)
			v1 = 2.*ran1(idum)-1.
			v2 = 2.*ran1(idum)-1.
			r = v1*v1 + v2*v2
		end do
		fac = sqrt(-2.*log(r)/r)
		gset = v1*fac
		iset = 1
		output = v2*fac
	else
		iset = 0
		output = gset
	end if
	return(output)
end

undef("weibul_params")
function weibul_params(skew[*][*]:float)
local sk,inva,Avec,Bvec,dims,ret_arr,skw,A
begin
	sk = (/-1.000,-0.971,-0.917,-0.867,-0.638,-0.254,0.069,0.359,0.631,0.896,1.160,1.430,1.708, \
			2.000,2.309,2.640,2.996,3.382,3.802,4.262,4.767,5.323,5.938,6.619,7.374,8.214/)
	inva = (/0.020,0.030,0.040,0.050,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800, \
			0.900,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.100,2.200/)
	Avec = (/0.446,0.444,0.442,0.439,0.425,0.389,0.346,0.297,0.246,0.193,0.142,0.092,0.044,0.000, \
			-0.040,-0.077,-0.109,-0.136,-0.160,-0.180,-0.196,-0.208,-0.217,-0.224,-0.227,-0.229/)
	Bvec = (/40.005,26.987,20.481,16.576,8.737,4.755,3.370,2.634,2.159,1.815,1.549,1.334,1.154,  \
			1.000,0.867,0.752,0.652,0.563,0.486,0.418,0.359,0.308,0.263,0.224,0.190,0.161/)

	dims = dimsizes(skew)
	ret_arr = new((/3,dims(0),dims(1)/),float)
	skw = where(skew.gt.sk(25),sk(25),skew)
	skw = where(skew.lt.sk(0),sk(0),skew)
	do ndx = 0,24
		A = (sk(ndx+1)-skw)/(sk(ndx+1)-sk(ndx))
		ret_arr(0,:,:) = where(skw.le.sk(ndx+1).and.skw.ge.sk(ndx),1/(A*inva(ndx)+(1-A)*inva(ndx+1)),ret_arr(0,:,:))
		ret_arr(1,:,:) = where(skw.le.sk(ndx+1).and.skw.ge.sk(ndx),A*Avec(ndx)+(1-A)*Avec(ndx+1),ret_arr(1,:,:))
		ret_arr(2,:,:) = where(skw.le.sk(ndx+1).and.skw.ge.sk(ndx),A*Bvec(ndx)+(1-A)*Bvec(ndx+1),ret_arr(2,:,:))
	end do
	return(ret_arr)
end

undef("get_data_from_F_weibul")
function get_data_from_F_weibul(mean[*][*]:float,stdev[*][*]:float,skew[*][*]:float,F[*][*]:float)
local b,vals,bound
begin
	vals = weibul_params(skew)
	b = vals(1,:,:)*stdev + mean
	bound = b - vals(2,:,:)*stdev
	bound = where(bound.lt.0,0.,bound)
	return((-log(1-F)^(1/vals(0,:,:)))*(b-bound) + bound)
end

undef("get_F_from_data_weibul")
function get_F_from_data_weibul(mean[*][*]:float,stdev[*][*]:float,skew[*][*]:float,x[*][*]:float,TINY[1]:float)
local vals,i,j,k,b,bound,F,diff
begin
	vals = weibul_params(skew)
	i = vals(0,:,:)
	j = vals(1,:,:)
	k = vals(2,:,:)
	b = j*stdev + mean
	bound = b - k*stdev
	bound = where(bound.le.0.or.bound.gt.x,0.,bound)
	diff = b-bound
	diff = where(diff.eq.0,diff@_FillValue,diff)
	F = 1 - exp(-1* ( ((x-bound)/diff) ^ i ) )
	F = where(F.ge.1.,F - TINY,F)
	F = where(ismissing(F).or.F.eq."nan".or.F.lt.TINY,TINY,F)
	return(F)
end

undef("get_data_from_F_EVI")
function get_data_from_F_EVI(mean[*][*]:float,stdev[*][*]:float,F[*][*]:float)
local a,b,fill
begin
	fill = 1e-20
	std = where(stdev.eq.0,fill,stdev)
	b = 3.14159/(std*sqrt(6))
	a = mean - 0.5772/b
	return(a - (1/b)*(log(-log(F))))
end

undef("get_F_from_data_EVI")
function get_F_from_data_EVI(mean[*][*]:float,stdev[*][*]:float,x[*][*]:float,TINY[1]:float)
local fill,std,a,b,F
begin
	fill = 1e-20
	std = where(stdev.eq.0,fill,stdev)
	b = 3.14159/(std*sqrt(6))
	a = mean - 0.5772/b
	F = exp(-1*exp(-b*(x-a)))
	F = where(F.ge.1.,F - TINY,F)
	F = where(F.lt.TINY,TINY,F)
	return(F)
end

undef("get_data_from_F_normal")
function get_data_from_F_normal(mean[*][*]:float,stdev[*][*]:float,F[*][*]:float)
local z
begin
	z = ((F^0.135)-((1-F)^0.135))/0.1975
	return(z*stdev + mean)
end

undef("get_F_from_data_normal")
function get_F_from_data_normal(mean[*][*]:float,stdev[*][*]:float,x[*][*]:float,TINY[1]:float)
local fill,std,z,F,sig
begin
	fill = 1e-20
	std = where(stdev.eq.0,fill,stdev)
	z = (x-mean)/std
	sig = where(z.lt.0,-1,1)
	z = abs(z)
	z = where(z.eq.0.,fill,z)
	F = 1-0.5*exp(-1*((83*z+351)*z+562)/(165+703/z))
	F = where(F.ge.1.,F-TINY,F)
	F = where(F.lt.TINY,TINY,F)
	F = where(sig.eq.-1,1-F,F)
	return(F)
end
	
undef("lookup")
function lookup(q[*][*]:float,v1[*][*][*]:float,v2[*][*][*]:float,par[1]:integer, \
				lu[1]:integer,mean[*][*]:float,stdev[*][*]:float,skew[*][*]:float)
local val,A,QUAN,PRCP,TEMP,DAT,TINY,dims,v1diff,vmask,F_weibul,F_EVI,F_norm
begin
	PRCP = 1
	TEMP = 0
	DAT = 1
	QUAN = 0
	dims = dimsizes(v1)
	TINY = 1./(dims(0)+1.)
	val = q
	val = 0.
	A = q
	A = 0.

	do ndx = 0,dims(0)-2
		v1diff = v1(ndx+1,:,:)-v1(ndx,:,:)
		v1diff = where(v1diff.eq.0,(v1(ndx+1,:,:)-q),v1diff)
		v1diff = where(v1diff.eq.0,-9999,v1diff)
		A = (v1(ndx+1,:,:)-q)/v1diff
		vmask = where(val.eq.0,True,False)
		val = where(q.le.v1(ndx+1,:,:).and.q.ge.v1(ndx,:,:).and.vmask,A*v2(ndx,:,:) + (1-A)*v2(ndx+1,:,:),val)
	end do
	delete(A)

	if (num(q.lt.v1(0,:,:)).gt.0.or.num(q.gt.v1(dims(0)-1,:,:)).gt.0) then
		if (lu.eq.QUAN) then
			if (par.eq.PRCP) then
				F_weibul = get_F_from_data_weibul(mean,stdev,skew,q,TINY)
				F_EVI = get_F_from_data_EVI(mean,stdev,q,TINY)
				val = where(q.lt.v1(0,:,:),F_weibul,val)
				val = where(q.gt.v1(dims(0)-1,:,:),F_EVI,val)
			end if
			if (par.eq.TEMP) then
				F_norm = get_F_from_data_normal(mean,stdev,q,TINY)
				val = where(q.lt.v1(0,:,:).or.q.gt.v1(dims(0)-1,:,:),F_norm,val)
			end if
		end if
		
		if (lu.eq.DAT) then
			if (par.eq.PRCP) then
				F_weibul = get_data_from_F_weibul(mean,stdev,skew,q)
				F_EVI = get_data_from_F_EVI(mean,stdev,q)
				val = where(q.lt.v1(0,:,:),F_weibul,val)
				val = where(q.gt.v1(dims(0)-1,:,:),F_EVI,val)
			end if
			if (par.eq.TEMP) then
				F_norm = get_data_from_F_normal(mean,stdev,q)
				val = where(q.lt.v1(0,:,:).or.q.gt.v1(dims(0)-1,:,:),F_norm,val)
			end if
		end if
	end if
	
	return(val)
end
