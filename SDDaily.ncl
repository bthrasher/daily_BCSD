load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
fill = 1.0e20
setfileoption("nc","Format","NetCDF4Classic") 
setfileoption("nc","CompressionLevel",5) ; 0 through 9 possible 

; Read in low-res BC files
bcfileshist = systemfunc("ls "+dhistbc+"/*.nc")	; files w/ BC factors
bcfiles21c = systemfunc("ls "+d21bc+"/*.nc")	; files w/ BC factors
bcfiles = array_append_record(bcfileshist,bcfiles21c,0)
nbcfiles = dimsizes(bcfiles)
bcfilesall = addfiles(bcfiles,"r")
if (dec.lt.dimsizes(bcfileshist)) then
	fout = str_sub_str(bcfiles(dec),dhistbc,dhistout)
else
	fout = str_sub_str(bcfiles(dec),d21bc,d21out)
end if	
if (isfilepresent(fout)) then exit end if
bcdims = dimsizes(bcfilesall[0]->$var$(:,{min_lat:max_lat},{min_lon:max_lon}))
nbclats = bcdims(1)
nbclons = bcdims(2)
delete(bcdims)
bclats = bcfilesall[0]->lat
bclons = bcfilesall[0]->lon

; Read in high-res obs file
ls_obsfiles = systemfunc("ls "+dobs+"/"+var+"*.nc | grep -v aggr")
obsfile = addfile(ls_obsfiles(0),"r")
delete(ls_obsfiles)
olats = obsfile->lat({min_lat:max_lat})
olons = obsfile->lon({min_lon:max_lon})
delete(obsfile)

; Create arrays
dectimes = bcfilesall[dec]->time
ts = dimsizes(dectimes)
decdates = cd_calendar(dectimes,-5)
if (decdates@calendar.eq."proleptic_gregorian") then decdates@calendar = "gregorian" end if
doy_dec = day_of_year(decdates(:,0),decdates(:,1),decdates(:,2)) - 1
oclimof = addfile(dobs+"/../"+var+"_obsclimo_"+ref_st+"-"+ref_end+".nc","r")
oclimo = oclimof->$var$(:,{min_lat:max_lat},{min_lon:max_lon})
opt = True
opt@critpc = 1
aggoclimo = area_hi2lores_Wrap(olons,olats,oclimo,True,1,bclons,bclats,opt)
if (var.eq."pr") then aggoclimo = where((aggoclimo*86400.).lt.0.25,0.,aggoclimo) end if

fmask = addfile(dobs+"/../../elevation_0.25.nc","r")
obsmask = fmask->elev(0,0,120:,:)
delete(fmask)

; Begin SD
print("Beginning SD")
do y = min(decdates(:,0)),max(decdates(:,0))
	print("year "+y)
	fstr = str_split(fout,"_")
	fouty = str_join(fstr(0:dimsizes(fstr)-2),"_") + "_" + y + ".nc"
	fouty = str_sub_str(fouty,"BC_"+var,var)
	delete(fstr)
	if (isfilepresent(fouty)) then 
		print("Output file "+fouty+" already exists")
		continue
	end if 
	ii = ind(decdates(:,0).eq.y)
	tcd = decdates(ii,:)
	t = bcfilesall[:]->time(ii)
	ts = dimsizes(t)
	doy_ii = day_of_year(tcd(:,0),tcd(:,1),tcd(:,2))

	factors = bcfilesall[dec]->$var$(ii,{min_lat:max_lat},{min_lon:max_lon})
	if (var.eq."pr") then factors = where((factors*86400.).lt.0.25,0.,factors) end if
	
	do d = 0,ts-1
		doy = doy_ii(d) - 1
		if (isStrSubset(var,"tas")) then
			factors(d,:,:) = factors(d,:,:) - aggoclimo(doy,:,:)
		else
			day_oclimo = aggoclimo(doy,:,:)
			day_oclimo = where(day_oclimo.eq.0,-2222,day_oclimo)			
			factors(d,:,:) = factors(d,:,:)/day_oclimo
			fmask = where(ismissing(day_oclimo),False,day_oclimo.eq.-2222)
			factors(d,:,:) = where(fmask,0.,factors(d,:,:))
		end if
	end do

	factors = linmsg(factors,-1)
	factors = linmsg_n(factors,-1,1)
	if (min_lon.eq.0.and.max_lon.eq.360) then
		factors_int = linint2_Wrap(bclons,bclats,factors,True,olons,olats,0)
	else
		factors_int = linint2_Wrap(bclons,bclats,factors,False,olons,olats,0)
	end if
	delete(factors)

	do d = 0,ts-1
		doy = doy_ii(d) - 1
		if (isStrSubset(var,"tas")) then
			factors_int(d,:,:) = factors_int(d,:,:) + oclimo(doy,:,:)
		else
			factors_int(d,:,:) = factors_int(d,:,:) * oclimo(doy,:,:)
		end if
		factors_int(d,:,:) = where(ismissing(obsmask),fill,factors_int(d,:,:))
	end do

	bcsd = addfile(fouty,"c")
	filedimdef(bcsd,"time",-1,True)
	bcsd->time = t
	bcsd->$var$ = factors_int
	delete([/fouty,ii,tcd,t,ts,doy_ii,factors_int/])
end do

end

