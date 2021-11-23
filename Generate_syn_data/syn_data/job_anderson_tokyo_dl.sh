suwaveform type=ricker1 fpeak=750 dt=0.0002  | sugain pbal=1 > ricker_wavelet.su
for loop in {1..500}
do
	radius=$(echo "${RANDOM}%10000*0.001+1"|bc)
  	./resonance_fre_yang ${radius}
  	label=$(echo "${radius}*10.00"|bc)
  	a2b < "${label}reflectivity_yang_tokyobaymm.txt" n1=5000 | suaddhead ns=5000 | sushw key=dt a=200 > "${label}reflectivity_yang_tokyobaymm.su"
  	suconv < "${label}reflectivity_yang_tokyobaymm.su" sufile=ricker_wavelet.su  | suwind tmin=0.0 tmax=0.99 | sustrip | suaddhead ns=5001 | sushw key=dt a=200 > "${label}waveform_yang_tokyobay.su"
  	sustrip < "${label}waveform_yang_tokyobay.su" | b2a n1=5001 > "${label}waveform_yang_tokyobay.txt"
  	rm "${label}waveform_yang_tokyobay.su"
done
