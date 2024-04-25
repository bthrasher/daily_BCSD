load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "BCFuncs.ncl"
load "get1Dindex2.ncl"

begin
pad = 4		; Defines detrend moving window, change as appropriate for your application
AVG_YRS = pad*2 + 1	; Number of years for moving average in BC detrend step
idum = 123
fill = 1.0e20
buff = 15	; Number of days +/- for CDF creation moving window

err = NhlGetErrorObjectId()
setvalues err
	"errLevel" : "Fatal"          ; only report Fatal errors
end setvalues

setfileoption("nc","Format","NetCDF4Classic") 
setfileoption("nc","CompressionLevel",5) ; 0 through 9 possible 

; Read in GCM files
files_hist = systemfunc("ls "+dhist+"/*.nc")
nfiles_hist = dimsizes(files_hist)
print("Found "+nfiles_hist+" 20C GCM files")

files21c = systemfunc("ls "+d21+"/*.nc")
nfiles21c = dimsizes(files21c)
print("Found "+nfiles21c+" 21C GCM files")

if (nfiles_hist.gt.1) then
	filesgcm = array_append_record(files_hist,files21c,0)
	nfiles = nfiles_hist + nfiles21c
else 
	filesgcm = files21c
	nfiles = nfiles21c
end if

; dec is the index of the file being adjusted in filesgcm
if (dec.lt.nfiles_hist) then
	fname = str_sub_str(filesgcm(dec),dhist,dhistout)
else
	fname = str_sub_str(filesgcm(dec),d21,d21out)
end if
if (isfilepresent(fname)) then 
	print(fname+" already exists. Exiting.")
	exit 
end if

fall = addfiles(filesgcm,"r")
gcmcal = str_lower(fall[0]->time@calendar)
print("GCM calendar = "+gcmcal)
if (isStrSubset(gcmcal,"360")) then
	ndays = 360
else
	ndays = 365
end if
startyrs = new((/nfiles/),integer)
endyrs = startyrs
do i = 0,nfiles-1
	t = fall[i]->time
	datest = cd_calendar(t(0),-5)
	dateend = cd_calendar(t(dimsizes(t)-1),-5)
	startyrs(i) = datest(0,0)
	endyrs(i) = dateend(0,0)
	delete(t)
	delete(datest)
	delete(dateend)
end do

fref_st = ind(startyrs.le.ref_st.and.endyrs.ge.ref_st)
fref_end = ind(startyrs.le.ref_end.and.endyrs.ge.ref_end)

; Calc 21C trend in temp data
if (isStrSubset(var,"tas")) then
	print("Calculating monthly trend")
	trend_times = yyyymm_time(startyrs(0),endyrs(nfiles-1),"integer")
	print("trend_times created")
	temp_mo_dat = new((/dimsizes(trend_times),dimsizes(fall[0]->lat({min_lat:max_lat})),dimsizes(fall[0]->lon({min_lon:max_lon}))/),float,fill)
	do f = 0,nfiles-1
		moAvg_temp = calculate_monthly_values(fall[f]->$var$(:,{min_lat:max_lat},{min_lon:max_lon}), "avg", 0, False)
		ii = ind(.not.ismissing(moAvg_temp(:,dimsizes(fall[0]->lat({min_lat:max_lat}))/2,dimsizes(fall[0]->lon({min_lon:max_lon}))/2)))
		moAvg = moAvg_temp(ii,:,:)
		delete(moAvg_temp)
		delete(ii)
		dims = dimsizes(moAvg)
		indst = ind(trend_times.eq.cd_calendar(fall[f]->time(0),-1))
		temp_mo_dat(indst:indst+dims(0)-1,:,:) = moAvg
		delete(moAvg)
		delete(dims)
	end do
	print("monthly avgs calculated")
	mon_avg_hist = clmMonTLL(temp_mo_dat(ind(trend_times.eq.(ref_st*100)+1):ind(trend_times.eq.(ref_end*100)+12),:,:))
	anom = calcMonAnomTLL(temp_mo_dat,mon_avg_hist)
	trend = new((dimsizes(temp_mo_dat)),float,fill)
	trend = 0.

	do m = 0,11
		all_mo = anom(m::12,:,:)
		temp = runave_n(all_mo,AVG_YRS,0,0)
		dims = dimsizes(temp)
		do p = 0,pad-1
			temp(dims(0)-pad+p,:,:) = temp(dims(0)-pad-1,:,:)
		end do
		trend(m::12,:,:) = temp
		delete(all_mo)
		delete(temp)
	end do
	
	trend(0:ind(trend_times.eq.(ref_end*100)+12),:,:) = 0.
	trend!0 = "time"
	trend!1 = "lat"
	trend!2 = "lon"
	
end if

; Read in aggregated obs files
ls_obsfiles = systemfunc("ls "+dobs+"/"+var+"_*.nc")
ii_stf = 0
ii_endf = dimsizes(ls_obsfiles)-1
do i = 0,dimsizes(ls_obsfiles)-1
	if (isStrSubset(ls_obsfiles(i),ref_st)) then ii_stf = i end if
	if (isStrSubset(ls_obsfiles(i),ref_end)) then ii_endf = i end if
end do
obsfiles = addfiles(ls_obsfiles(ii_stf:ii_endf),"r")
lats = obsfiles[0]->lat({min_lat:max_lat})
lons = obsfiles[0]->lon({min_lon:max_lon})
nlats = dimsizes(lats)
nlons = dimsizes(lons)
timesobs = obsfiles[:]->time
datesobs = cd_calendar(timesobs,-5)
delete(timesobs)
doyobs = day_of_year(datesobs(:,0),datesobs(:,1),datesobs(:,2)) - 1
delete(datesobs)

; Create time arrays
dectimes = fall[dec]->time
ts = dimsizes(dectimes)
decdates = cd_calendar(dectimes,-5)
if (decdates@calendar.eq."proleptic_gregorian") then decdates@calendar = "gregorian" end if
dec_yyyymm = cd_calendar(dectimes,-1)
dates_hist_temp = cd_calendar(fall[fref_st:fref_end]->time,-5)
ii_dates = ind(dates_hist_temp(:,0).ge.ref_st.and.dates_hist_temp(:,0).le.ref_end)
dates_hist = dates_hist_temp(ii_dates,:)
if (dates_hist@calendar.eq."proleptic_gregorian") then dates_hist@calendar = "gregorian" end if
delete(dates_hist_temp)
ts_hist = dimsizes(dates_hist(:,0))
doy_hist = day_of_year(dates_hist(:,0),dates_hist(:,1),dates_hist(:,2)) - 1
doy_dec = day_of_year(decdates(:,0),decdates(:,1),decdates(:,2)) - 1
delete(dates_hist)

; Get daily GCM data from decadal file
decdata = fall[dec]->$var$(:,{min_lat:max_lat},{min_lon:max_lon})
if (typeof(decdata).eq."double") then decdata = doubletofloat(decdata) end if
if (isStrSubset(var,"tas")) then
	print("Removing monthly trend from daily data")
	do d = 0,ts-1
		itrend = ind(trend_times.eq.dec_yyyymm(d))
		decdata(d,:,:) = decdata(d,:,:) - trend(itrend,:,:)
	end do
end if

obsmask = where(ismissing(obsfiles[0]->$var$(0,{min_lat:max_lat},{min_lon:max_lon})),False,True)
bc_data_all = new((/ts,nlats,nlons/),float,fill)

print("Correcting bias in dec "+dec)

do d = 0,ts-1
	doy = doy_dec(d)
	
	print("Correcting bias in dec "+dec+", year "+decdates(d,0)+", month "+decdates(d,1)+", day "+decdates(d,2)+ ", doy = "+doy+", ts = "+d)
	
	if (doy.gt.ndays) then doy = ndays end if
	doybuff = ispan(doy-buff,doy+buff,1)
	doybuff_hist = doybuff
	doybuff_hist = where(doybuff_hist.lt.0,doybuff_hist+ndays,doybuff_hist)
	doybuff_hist = where(doybuff_hist.ge.ndays,doybuff_hist-ndays,doybuff_hist)
	doybuffobs = doybuff
	doybuffobs = where(doybuffobs.lt.0,doybuffobs+ndays,doybuffobs)
	doybuffobs = where(doybuffobs.ge.365,doybuffobs-ndays,doybuffobs)
	
	ii = get1Dindex2(doy_hist,doybuff_hist)
	qsort(ii)
	temp_hist = fall[fref_st:fref_end]->$var$(ii_dates,{min_lat:max_lat},{min_lon:max_lon})
	sort_hist = temp_hist(ii,{min_lat:max_lat},{min_lon:max_lon})
	delete(temp_hist)
	delete(ii)

	ii = get1Dindex2(doyobs,doybuffobs)
	qsort(ii)
	sortobs = obsfiles[:]->$var$(ii,{min_lat:max_lat},{min_lon:max_lon})
	delete(ii)

	dims_hist = dimsizes(sort_hist)
	dimsobs = dimsizes(sortobs)
	do i = 0,dims_hist(0)-1
		sort_hist(i,:,:) = where(obsmask,sort_hist(i,:,:),fill)
	end do
	quant_hist = new((/dims_hist(0),nlats,nlons/),float,fill)	
	do i = 0,dims_hist(0)-1
		quant_hist(i,:,:) = (i+1.)/(dims_hist(0)+1.)
	end do
	quantobs = new((/dimsobs(0),nlats,nlons/),float,fill)	
	do i = 0,dimsobs(0)-1
		quantobs(i,:,:) = (i+1.)/(dimsobs(0)+1.)
	end do
	
	day_data = decdata(d,:,:)	
	day_data = where(obsmask,day_data,fill)
	
	print("Sorting 20C and obs")
	ip = dim_pqsort_n(sortobs,2,0)
	delete(ip)
	ip = dim_pqsort_n(sort_hist,2,0)
	delete(ip)
	
;	Perform stats
	stat_hist = dim_stat4_n(sort_hist,0)
	mean_hist = stat_hist(0,:,:)
	var_hist = stat_hist(1,:,:)
	stdev_hist = sqrt(var_hist)
	skew_hist = stat_hist(2,:,:)
	skew_hist = where(var_hist.eq.0,0.,skew_hist)
	delete(stat_hist)
	delete(var_hist)
	statobs = dim_stat4_n(sortobs,0)
	meanobs = statobs(0,:,:)
	varobs = statobs(1,:,:)
	stdevobs = sqrt(varobs)
	skewobs = statobs(2,:,:)
	skewobs = where(varobs.eq.0,0.,skewobs)
	delete(statobs)
	delete(varobs)
	
	print("starting lookup")
	if (var.eq."pr") then
		q = lookup(day_data,sort_hist,quant_hist,1,0,mean_hist,stdev_hist,skew_hist)
		bc_data = lookup(q,quantobs,sortobs,1,1,meanobs,stdevobs,skewobs)
		bc_data = where(bc_data.le.0.,0.,bc_data)
	else
		q = lookup(day_data,sort_hist,quant_hist,0,0,mean_hist,stdev_hist,skew_hist)
		bc_data = lookup(q,quantobs,sortobs,0,1,meanobs,stdevobs,skewobs)
	end if
	
	bc_data_all(d,:,:) = bc_data
	
	delete(sort_hist)
	delete(sortobs)
	delete(q)
	delete(bc_data)
	delete(quant_hist)
	delete(quantobs)
		
end do	; End day loop

	
if (isStrSubset(var,"tas")) then
	print("Returning monthly trend to daily data")
	do d = 0,ts-1
		itrend = ind(trend_times.eq.dec_yyyymm(d))
		bc_data_all(d,:,:) = bc_data_all(d,:,:) + trend(itrend,:,:)
	end do
end if
	
if (isatt(dectimes,"bounds")) then delete(dectimes@bounds) end if
if (isatt(lats,"bounds")) then delete(lats@bounds) end if
if (isatt(lons,"bounds")) then delete(lons@bounds) end if

bc_data_all!0 = "time"
bc_data_all!1 = "lat"
bc_data_all!2 = "lon"
bc_data_all&time = dectimes
bc_data_all&lat = lats
bc_data_all&lon = lons

print("Creating file "+fname)
fout = addfile(fname,"c")
filedimdef(fout,"time",-1,True)
fout->$var$ = bc_data_all
fout->time = dectimes
fout->lat = lats
fout->lon = lons

end


