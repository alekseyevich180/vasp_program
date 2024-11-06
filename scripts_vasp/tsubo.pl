#!/usr/bin/perl
#use strict;
my $version="Tsuboguruma by Hiroaki Nissho Ver 2681/2/21\n";
#my $phonopy_tolerance=1E-2;
my $phonopy_tolerance=1E-4;
#POSCAR utility Tsuboguruma 
#Original version by Hiroaki Nissho
#NO GUARANTEE/WARRANTY FOR USE OF THIS CODE
#Support by Hiroaki Nissho provided only in Japanese -- no exceptions

#Saisin_ban ha http://tsuboguruma.seesaa.net/ wo sansyou

#kore wa tsukawanai
my @arg=@ARGV;
my @arg2=@ARGV;
shift @arg2;
my ($topline,$scale,@latvec,@namespecies,@numspecies,$num_atoms,$dc,@x,@y,@z,@w);

# A1) Help
&help if ($ARGV[0] eq "-h");
&help if ($ARGV[0] eq "-help");
&help if ($ARGV[0] eq "--help");
&version if ($ARGV[0] eq "--ver");
&version if ($ARGV[0] eq "--version");

# A2) No POSCAR, sonota
if ($ARGV[0] eq "--clat"){
	&clat(@arg2[0..5]);
}

if (($ARGV[0] eq "--symop_info") || ($ARGV[0] eq "--symmetry_operation_info")) {
	my @a=&get_isometry_type(@arg2[0..11]);
	print ("$a[0] [ @a[1..3] ] ( @a[4..12] ) ( @a[13..15] )\n");
	exit;
}



#Read POSCAR
&readPOSCAR;

#dummy: for debugging
if ($ARGV[0] eq "--dummy") {
	exit;
}



# B) POSCAR info

if ($ARGV[0] eq "--data") {
	&data("normal");
}

if ($ARGV[0] eq "--data_full") {
	&data("full");
}

if (($ARGV[0] eq "--dataniggli") || ($ARGV[0] eq "--data_niggli")){
	@latvec=&get_niggli_latvec(@latvec);
	&data;
}

if (($ARGV[0] eq "--natoms") || ($ARGV[0] eq "--numatoms") || ($ARGV[0] eq "--num_atoms") || ($ARGV[0] eq "--natom") || ($ARGV[0] eq "--numatom") || ($ARGV[0] eq "--num_atom")){
	print ("$num_atoms\n");
	exit;
};

if (($ARGV[0] eq "--nspecies") || ($ARGV[0] eq "--numspecies") || ($ARGV[0] eq "--num_species")){
	my $a=@numspecies;
	print ("$a\n");
	exit;
}

if ($ARGV[0] eq "--stoichiometry"){
	my @a=&minimum_array(@numspecies);
	print ("@a\n");
	exit;
}

if ($ARGV[0] eq "--which_atom_number"){
	my $a=&which_atom_number(@arg2[0..2]);
	print ("$a\n");
	exit;
}


if (($ARGV[0] eq "--dist")||($ARGV[0] eq "--distance")){
	&distance(1,"default",$arg2[0]);
}

if (($ARGV[0] eq "--dist_large")||($ARGV[0] eq "--distance_large")){
	&distance(3,"default",$arg2[0]);
}

if (($ARGV[0] eq "--dist_min")||($ARGV[0] eq "--distance_min")){
	&distance(1,"min");
}

if (($ARGV[0] eq "--dist_atom_list")||($ARGV[0] eq "--distance_atom_list")){
	shift @arg2;
	my @atoms;
	my $a;
	while ($a=shift (@arg2)){
                my @b=split (/\.\./, $a);
                $b[1]=$b[0] if ($b[1] eq "");
                @atoms = (@atoms, $b[0] .. $b[1]);
	}
	&distance(1,"list",$ARGV[1], @atoms);
}

if ($ARGV[0] eq "--solid_angle"){
	&distance_sphere($ARGV[1]);
}

if (($ARGV[0] eq "--coord") || ($ARGV[0] eq "--coordination")){
	&distance(1,"coordination",@arg2[0..1]);
}

if (($ARGV[0] eq "--coordination_large") || ($ARGV[0] eq "--coord_large")){
	&distance(2,"coordination",@arg2[0..1]);
}

if (($ARGV[0] eq "--maxgap") || ($ARGV[0] eq "--max_gap")) {
#z_max_gap
	&c2d if ($dc eq "C");
	my @copy;
	foreach my $i (@z){
		$i=$i-&floor($i);
		push (@copy, ($i*0.5, $i*0.5+0.5));
	}
	my @gapinfo=&max_gap(@copy);
	my $gapsize=$gapinfo[1]*2;

#vol/|axb|
	my @ab=&cross_vv(@latvec[0..5]);
	my $vol=&product_vv(@ab,@latvec[6..8])*$scale;
	my $cnorm=$vol/&norm(@ab)*$scale;
	my $gapsizecart=$gapsize*$cnorm;
	print("gap_size $gapsize $gapsizecart \n");
	exit;
}

if (($ARGV[0] eq "--volperatom") || ($ARGV[0] eq "--vol_per_atom")){
	my $density=&det3(@latvec)*$scale*$scale*$scale/$num_atoms;
	print ("$density \n");
	exit;
}


if (($ARGV[0] eq "--nameatom")||($ARGV[0] eq "--name_atom")){
#Genso mei to bango wo hyoji
	my @a=&nameatom($arg2[0]);
	print ("@a \n");
	exit;
}

if (($ARGV[0] eq "--kpoints") || ($ARGV[0] eq "--kpoints_slab")) {
	my @info=&kpointscheck;
	for (my $i=0; $i<=2; $i++){
		my $a=int($info[$i]/$ARGV[1]+0.7);
		$info[$i]=&max($a,1);
	}
	$info[2]=1 if ($ARGV[0] eq "--kpoints_slab");
	&printkpoints(@info);
}

if (($ARGV[0] eq "--kpoints_even") || ($ARGV[0] eq "--kpoints_even_slab")) {
	my @info=&kpointscheck;
	for (my $i=0; $i<=2; $i++){
		my $a=2*int($info[$i]/$arg2[0]/2+0.7);
		$info[$i]=&max($a,2);
	}
	$info[2]=1 if ($ARGV[0] eq "--kpoints_even_slab");
	&printkpoints(@info);
}

sub kpointscheck{
#gyakukousi a b c to centering wo kaesu
	die ("mesh kankaku $arg2[0] ga arimasen \n") if ($ARGV[1] eq "");
	my $center="G";
	if ($arg2[1] ne ""){
		if ($arg2[1]=~/^M/){
			$center="M";
		} elsif ($arg2[1]=~/^G/){
		} else {
			die ("center type $arg2[1] ga M ka G de hajimarimasen");
		}
	}
	my @recilatvec=&reci_latvec(@latvec);
	my @abcreci=&get_abc(@recilatvec);
	return (@abcreci[0..2],$center);
}

sub printkpoints{
#print kpoints
	print ("Automatic mesh \n");
	print ("0 \n");
	print ("$_[3] \n");
	print ("@_[0..2] \n");
	print ("0 0 0 \n");
	exit;
}

if (($ARGV[0] eq "--slabarea")||($ARGV[0] eq "--slab_area")){
	my $area=2*&norm(&cross_vv(@latvec[0..5]));
	print ("$area \n");
	exit;
}

if (($ARGV[0] eq "--average_atom")||($ARGV[0] eq "--average_site")){
	my ($x_ave,$y_ave,$z_ave,$count);
	while ($arg2[0] ne ""){
		$count++;
		$x_ave+=$x[$arg2[0]-1];
		$y_ave+=$y[$arg2[0]-1];
		$z_ave+=$z[$arg2[0]-1];
		shift @arg2;
	}
	$x_ave/=$count;
	$y_ave/=$count;
	$z_ave/=$count;
	print ("$x_ave $y_ave $z_ave \n");
	exit;
}


if (($ARGV[0] eq "--sgnumber")||($ARGV[0] eq "--sg_number")){
	my @a=&get_bposcar;
	print ("$a[0]\n");
	exit;
}

if (($ARGV[0] eq "--sgsymbol")||($ARGV[0] eq "--sg_symbol")){
	my @a=&get_bposcar;
	print ("$a[1]\n");
	exit;
}

if ($ARGV[0] eq "--2Dsgnumber"){
	my @a=&twoDspacegroup($arg2[0]);
	print ("$a[0]\n");
	exit;
}

if ($ARGV[0] eq "--2Dsgsymbol"){
	my @a=&twoDspacegroup($arg2[0]);
	print ("$a[1]\n");
	exit;
}

if ($ARGV[0] eq "--pgsymbol"){
	my @a=&get_bposcar;
	print ("$a[6]\n");
	exit;
}

if (($ARGV[0] eq "--crystalsystem") || ($ARGV[0] eq "--crystal_system")){
	my @a=&get_bposcar;
	print ("$a[3]\n");
	exit;
}

if ($ARGV[0] eq "--bravais"){
	my @a=&get_bposcar;
	print ("$a[4]\n");
	exit;
}

if ($ARGV[0] eq "--extended_bravais"){
	my @a=&get_bposcar;
	print ("$a[5]\n");
	exit;
}

if (($ARGV[0] eq "--sclatticetype") || ($ARGV[0] eq "--latticetype")){
	my $a=&get_standard_SC("conventional");
	print ("$a\n");
	exit;
}

#HKL hyoumen de mennai ni suichoku na nagasa $arg2[3] no vecto wo hyouji
if ($ARGV[0] eq "--normal_vector"){
	&get_hkl_cell("primitive",@arg2[0..2]);
	my @c=&cross_vv(@latvec[0..5]);
	@c=&product_vs(@c,$arg2[3]/&norm(@c));
	&printx($latvec[0]);
	&printx($latvec[1]);
	&printx($latvec[2]);
	print ("\n");
	&printx($latvec[3]);
	&printx($latvec[4]);
	&printx($latvec[5]);
	print ("\n");
	&printx($c[0]);
	&printx($c[1]);
	&printx($c[2]);
	print ("\n");
	exit;
}


if ($ARGV[0] eq "--isometry"){
	my @all_isometry=&get_isometry;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		print ("@all_isometry[$i..($i+11)]");
		if ($ARGV[1] eq "info"){
			my @a=&get_isometry_type(@all_isometry[$i..($i+11)]);
			print (" type $a[0] [ @a[1..3] ] ( @a[4..12] ) ( @a[13..15] )");
		}
		print ("\n");
	}
	exit;
}



if (($ARGV[0] eq "--isometry_nonpolar") || ($ARGV[0] eq "--isometry_nonpolar_image") || ($ARGV[0] eq "--isometry_nonpolar_square")){
#primitive wo tsukuru
	my @a=&get_isometry;
	my $count=0;
	for (my $i=0; $i<$#a; $i+=12){
		if (&norm(&diff_vv(@a[($i+6)..($i+8)],0,0,-1)) == 0){
			$count++;
			if ($ARGV[0] eq "--isometry_nonpolar"){
				print ("@a[$i..($i+11)]");
				if ($ARGV[1] eq "info"){
					my @b=&get_isometry_type(@a[$i..($i+11)]);
					print (" type $b[0] [ @b[1..3] ] ( @b[4..12] ) ( @b[13..15] )");
				}
				print ("\n");
#--isometry_nonpolar_image
#site bangou ka xyz wo hikisuu to site watasu
#kaerichi ha site bangou : site bangou , xyz : xyz
			} elsif ($ARGV[0] eq "--isometry_nonpolar_image"){
				if ($dc eq "C"){
					$dc="D";
					&c2d;
				}
				my @image;
				if ($ARGV[2] eq ""){
					@image=&isometry_image(@a[$i..($i+11)],$ARGV[1]);
				} else {
					@image=&isometry_image(@a[$i..($i+11)],@ARGV[1..3]);
				}
				print (" @image for @a[$i..($i+11)]\n");
#--isometry_nonpolar_square
			} else {
				my @square=&product_symop(@a[$i..($i+11)],@a[$i..($i+11)]);
				my $dist=&norm(&diff_vv(@square[0..8],1,0,0,0,1,0,0,0,1));
				if ($dist < 0.001){
					my @b=&isometry_real(@a[$i..($i+11)]);
					print ("@b translation of square: @square[9..11] \n");
				}
			}
		} 
	}
	print ("n/a \n") if ($count == 0);
	exit;
}


if ($ARGV[0] eq "--surface_isometry"){
#primitive wo tsukuru
	&get_hkl_cell("primitive",@arg2[0..3]);
	my @a=&get_surface_isometry;
	if ($#a > 0){
		for (my $i=0; $i<$#a; $i+=12){
			print ("@a[$i..($i+11)] \n");
		}
	} else {
		print ("n/a \n");
	}
	exit;
}


if ($ARGV[0] eq "--unique_nonpolar"){
#crystallographic conventional ni suru
#cap: menseki ga saisyou slab no $cap bai wa hyouji sinai
	my $cap=9999;
	$cap=$arg2[3] if ($arg2[3] ne "");
	my @a=&get_bposcar;
	my @data;
	my ($hmax, $kmax, $lmax)=(4,4,4);
	$hmax=$arg2[0] if ($arg2[0] ne "");
	$kmax=$arg2[1] if ($arg2[1] ne "");
	$lmax=$arg2[2] if ($arg2[2] ne "");
	for (my $h=0; $h<=$hmax; $h++){
		for (my $k=-$kmax; $k<=$kmax; $k++){
			for (my $l=-$lmax; $l<=$lmax; $l++){
				push @data, ($h, $k, $l) if (&is_unique_nonpolar($a[0],$h,$k,$l));
			}
		}
	}
	die ("nonpolar suface ga nai: space group $a[0] \n") if ($data[0] eq "");
	&nonpolar_area_sort($cap,@data);
	exit;
}

if ($ARGV[0] eq "--to_unique_nonpolar"){
#unique nonpolar wo kaesu
	my @a=&get_bposcar;
	my @temp=&to_unique_nonpolar($a[0],@arg2[0..2]);
	print ("@temp");
	exit;
}

if ($ARGV[0] eq "--isometry_nonpolar_image_file"){
#primitive wo tsukuru
	my @a=&get_isometry;
	my $count=0;
	my $tolerance=1E-4;
	if (! -e $arg2[0]){
		die ("--isometry_nonpolar_image_file: file $arg2[0] ga sonzai shinai\n");
	} else {
		open A, $arg2[0];
		my @decimal_arg=&diff_vv(@arg2[10..12],&floor_array(@arg2[10..12]));
#true isometry check
		for (my $i=0; $i<$#a; $i+=12){
			my @symop=@a[$i..($i+11)];
#tsubo.pl no vector part ha 0 to 1 no aida
			if (&norm(&diff_vv(@a[$i..($i+11)],@arg2[1..9],@decimal_arg)) < $tolerance){
#true isometry wo nyuusyu
				my @isometry_use=&isometry_real(@a[$i..($i+11)]);
#tadasii vector part
				@isometry_use[9..11]=&sum_vv(@isometry_use[9..11],&floor_array(@arg2[10..12]));
#matrix+vector
				while (my $b=<A>){
					my @site=&splitall($b);
					my @c=&product_mv(@isometry_use[0..8],@site);
					@c=&sum_vv(@isometry_use[9..11],@c);
					print ("@c \n");
				}
				exit;
			}
		}
	}
	print ("--isometry_nonpolar_image_file: isometry @arg2[1..12] ga sonzai shinai\n");
	exit;
}



#hyoumen no kotonaru termination no polarity no data wo kaesu: osusume
if ($ARGV[0] eq "--termination_polarity"){
	my @a=&get_termination_polarity(@arg2);
	my $t=`grep -c Nonpolar temp_tsubo_polarity_data`*1;
	if ($t != 0){
		system ("grep -v Non-stoichiometry temp_tsubo_polarity_data | grep -v Semipolar");
	} else {
		system ("grep -v Non-stoichiometry temp_tsubo_polarity_data");
	}
	system ("rm temp_tsubo_polarity_data");
	exit;
}

#hyoumen no kotonaru termination no polarity no data wo kaesu
if ($ARGV[0] eq "--termination_polarity_full"){
	&get_termination_polarity(@arg2);
	system ("cat temp_tsubo_polarity_data");
	system ("rm temp_tsubo_polarity_data");
	exit;
}


#hyoumen no kotonaru termination no polarity no data wo kaesu:

if ($ARGV[0] eq "--termination_polarity_list"){
	die ("file $arg2[0] ga sonzai sinai\n") if (! -e $arg2[0]);
	open LIST, "$arg2[0]";
	print ("H K L Polarity \n");
	while (my $t=<LIST>){
		my @a=&splitall($t);
		my @b=&get_termination_polarity(@a[0..2]);
		print ("@a[0..2] $b[0] @a[3..$#a+1] $b[5] \n");
	}
	close LIST;
	system ("rm temp_tsubo_polarity_data");
	exit;
}


#inplane_lattice_vector wo sort site kaesu
if ($ARGV[0] eq "--inplane_lattice_vector"){
	my $a=&get_inplane_lattice_vector(@arg2[0..3]);
	exit;
}


#terrace, vicinal vector wo sort site kaesu
if (($ARGV[0] eq "--terrace_vicinal_vector") || ($ARGV[0] eq "--vicinal_terrace_vector")){
	my $a=&get_terrace_vicinal_vector(@arg2[0..6]);
	exit;
}

#orientation no soutai haichi wo siraberu
if (($ARGV[0] eq "--orientation_relative_all") || ($ARGV[0] eq "--relative_orientation_all")){
	die ("orientation_relative_all: hikisuu ga 6 nai") if ($arg2[5] eq "");
	&facet_orientation(@arg2[0..5]);
	exit;
}

#orientation no soutai haichi wo siraberu: facet only
if (($ARGV[0] eq "--orientation_relative_facet") || ($ARGV[0] eq "--relative_orientation_facet")){
	die ("orientation_relative_facet: hikisuu ga 6 nai") if ($arg2[5] eq "");
	my @a=&get_bposcar;
	&facet_orientation(@arg2[0..5],$a[0]);
	exit;
}


#facet reconstruction surface energy wo siraberu
if (($ARGV[0] eq "--facet_surface_energy") || ($ARGV[0] eq "--facet_energy") || ($ARGV[0] eq "--surface_energy_facet")){
	die ("facet_surface_energy: hikisuu ga 7 nai") if ($arg2[6] eq "");
	die ("facet_surface_energy: library file $arg2[6] ga nai") if (! -e $arg2[6]);
	my @a=&get_bposcar;
	&facet_orientation(@arg2[0..5],$a[0], $arg2[6]);
	exit;
}

if (($ARGV[0] eq "--stepenergy") || ($ARGV[0] eq "--step_energy")){
#step energy density beta wo kaku
#nyuuryokku:  vicinal surface energy, terrace surface energy,  gensi bangou
	my (@newx,@newy,@newz);
	my $vicE= shift @arg2;
	my $terE= shift @arg2;
	foreach my $i (@arg2){
		push (@xnew,$x[$i-1]);
		push (@ynew,$y[$i-1]);
		push (@znew,$z[$i-1]);
	}
#z de sort
	my @order_key=&order_key_number(@znew);
	@x=@xnew[@order_key];
	@y=@ynew[@order_key];
	@z=@znew[@order_key];
	$num_atoms=$#x+1;
#latvec wo keru
	my @abc=&get_abc(@latvec);
	my @clat=&clat(@abc,"s");
	@latvec=&product_mm($clat[0],0,0,@clat[1..2],0,@clat[3..5],&sind($abc[5]),&cosd($abc[5]),0,-&cosd($abc[5]),&sind($abc[5]),0,0,0,1);

#cartesian
	&d2c;
#x wo seiri
	foreach (my $i=0; $i<$num_atoms;$i++){
		my $xpos=($x[$i]-$x[$num_atoms-1])/$latvec[0]+0.05;
		$x[$i]-=$latvec[0]*&floor($xpos);
	}
	#least squares
	my @lsq=&leastsquares(@x,@z);
	$angle=atan2(-$lsq[0],1);

	my $beta=($vicE-$terE* cos($angle))/sin($angle);
	my $height=$latvec[0]* sin($angle);
	my $beta_height=$beta*$height;
	print ("vicinal_Esurf $vicE terrace_Esurf $terE beta_normalized $beta beta_times_h $beta_height height $height theta(rad) $angle \n");
	exit;
}


if ($ARGV[0] eq "--equivalent_site") {
	die ("gensi $arg2[0] ga sonzai sinai \n") if ($arg2[0] > $num_atoms);
#to d
	my @isometry_print;
	my @equivalent_print;
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
	my @latvec_orig=@latvec;
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @numspecies_orig=@numspecies;
	my $num_atoms_orig=$num_atoms;
	my @a=&get_bposcar;
	print ("Kuukan gun @a[0..1] \n");
	@latvec=@latvec_orig;
	@x=@x_orig;
	@y=@y_orig;
	@z=@z_orig;
	@numspecies=@numspecies_orig;
	$num_atoms=$num_atoms_orig;
#isometry list
	my @all_isometry=&get_isometry;
#siraberu gensi
	my @atom=($x[$arg2[0]-1],$y[$arg2[0]-1],$z[$arg2[0]-1]);
#touka gensi list
	my @equivalent_atom_list;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		push (@equivalent_atom_list, &near_atom(&apply_symop(@atom,@all_isometry[$i..($i+11)])));
	}
	@a=&unique_array(@equivalent_atom_list);
	my @order_key=&order_key_number(@a);
	@equivalent_atom_list=@a[@order_key];
	print ("Touka na gensi: @equivalent_atom_list \n");
	for (my $i=0;$i<$#all_isometry;$i+=12){
#taisyou sousa I wo jyogai
		if (&norm(&diff_vv(@all_isometry[$i..($i+11)],1,0,0,0,1,0,0,0,1,0,0,0)) > 0.001){
			my @atom_set;
			my @living_atom_list=@equivalent_atom_list;
			my @isometry_check=@all_isometry[$i..($i+11)];
			my @symop_type=&get_isometry_type(@isometry_check);
			my $isometry_data=join (" ",@isometry_check);
			push (@isometry_print, "Type ".$symop_type[0]." : ".$isometry_data." :");
#set wo tsukuru
#ikinokori wo kesiteiku
			while (my $atom_num=shift (@living_atom_list)){
				my @atom=($x[$atom_num-1],$y[$atom_num-1],$z[$atom_num-1]);

				my $flag=1;
				my @isometry=@isometry_check;
				my @same_atom_list = ($atom_num);
#onaji isometry de tsunagaru gensi wo dasu
				while ($flag == 1){
					my $j=&near_atom(&apply_symop(@atom,@isometry));
					push (@same_atom_list,$j);
					@isometry=&product_symop(@isometry_check,@isometry);
					@isometry[9..11]=&diff_vv(@isometry[9..11],&round_array(@isometry[9..11]));
					$flag=0 if (&norm(&diff_vv(@isometry,1,0,0,0,1,0,0,0,1,0,0,0)) < 0.001);
				}
#tsunagaru gensi wo _ de tsunageru
				@a=&unique_array(@same_atom_list);
				my @order_key=&order_key_number(@a);
				@same_atom_list=@a[@order_key];
				push (@atom_set, join ("_",@same_atom_list));
#tsukatta gensi wo kesu
				foreach my $k (@same_atom_list) {
					@living_atom_list = grep $_ ne $k, @living_atom_list;
				}
			}
			push (@equivalent_print, join (" ",@atom_set));
		}
	}
	if (($arg2[1] ne "no_sort") && ($arg2[1] ne "nosort")){
		my @order_key=&order_key_cmp(@equivalent_print);
		@equivalent_print=@equivalent_print[@order_key];
		@isometry_print=@isometry_print[@order_key];
	}
	for (my $i=0; $i<=$#isometry_print; $i++){
		print ("$isometry_print[$i] $equivalent_print[$i] \n");
	}
	exit;
}




if ($ARGV[0] eq "--2ddiffraction"){
#check
	my @hkl_orig=@arg2[0..2];
	my $count=0;
	for (my $i=0;$i<=2;$i++){
		die ("seisuu igai ga @hkl_orig ni aru\n") if ($hkl_orig[$i] !~ /^[+-]?[0-9]+$/);
		$count+=($i+1)*($i+1) if ($hkl_orig[$i] == 0);
	}
	die ("@hkl_orig ga zenbu 0 \n") if ($count==14);
#2d lattice vector
	my @hkl1=(-$hkl_orig[1],$hkl_orig[0],0);
	my @hkl2=(-$hkl_orig[2],0,$hkl_orig[0]);
	my @hkl3=(0,-$hkl_orig[2],$hkl_orig[1]);
	@hkl1=&gcm_direction(@hkl1);
	@hkl2=&gcm_direction(@hkl2);
	@hkl3=&gcm_direction(@hkl3);
	@hkl1=@hkl3 if (&product_vv(@hkl1,@hkl1) == 0);
	@hkl2=@hkl3 if (&product_vv(@hkl2,@hkl2) == 0);
#011, 0-1-1 reigai
	my @deltahkl12=&diff_vv(@hkl1,@hkl2);
	@hkl2=@hkl3 if (&product_vv(@deltahkl12,@deltahkl12) == 0);
	@deltahkl12=&sum_vv(@hkl1,@hkl2);
	@hkl2=@hkl3 if (&product_vv(@deltahkl12,@deltahkl12) == 0);
#@hkl1,@hkl2 wo tsukau: migitekei ni suru
	@hkl1=&gcm_direction(@hkl1);
	@hkl2=&gcm_direction(@hkl2);
	my $s=&det3(@hkl_orig,@hkl1,@hkl2);
	@hkl2=&product_vs(@hkl2,-1) if ($s<0);
#Gaussian lattice reduction no jyunbi
#lattice vector: @hklreci1,2
#houi: @hkl1,2
	my @latvec_reci=&tmatrix3(&invmatrix3(@latvec));
	my @hklreci1=&product_vm(@hkl1,@latvec_reci);
	my @hklreci2=&product_vm(@hkl2,@latvec_reci);
	my @reduced=&gaussian_reduction_2d(@hklreci1,@hklreci2,@hkl1,@hkl2);
	@hklreci1=@reduced[0..2];
	@hklreci2=@reduced[3..5];
	@hkl1=@reduced[6..8];
	@hkl2=@reduced[9..11];
#Point 3: @hkl2-@hkl1
	my @hklreci3=&diff_vv(@hklreci2,@hklreci1);
	@hkl3=&diff_vv(@hkl2,@hkl1);
#Kyori
	my @distance=(1,&norm(@hklreci2)/&norm(@hklreci1),&norm(@hklreci3)/&norm(@hklreci1));
#Kakudo
	my @cosines=(1,&product_vv(@hklreci1,@hklreci2)/&norm(@hklreci1)/&norm(@hklreci2),&product_vv(@hklreci1,@hklreci3)/&norm(@hklreci1)/&norm(@hklreci3));
	my @angles=(0,&acos($cosines[1]),&acos($cosines[2]));
	my @x_coord=(1,$distance[1]*$cosines[1],$distance[2]*$cosines[2]);
	$cosines[1]=&precise($cosines[1],0);
	my @sines=(0,sqrt(1-$cosines[1]*$cosines[1]),sqrt(1-$cosines[2]*$cosines[2]));
	my @y_coord=(0,$distance[1]*$sines[1],$distance[2]*$sines[2]);
	print ("Diffraction spots\n");
	print ("Index Kyori Kakudo X-zahyo Y-zahyo\n");
	printf ("%3d %3d %3d ",@hkl1);
		printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[0], $angles[0], $x_coord[0], $y_coord[0]);
	printf ("%3d %3d %3d ",@hkl2);
	printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[1], $angles[1], $x_coord[1], $y_coord[1]);
	if ($cosines[1] != 0){
		printf ("%3d %3d %3d ",@hkl3);
		printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[2], $angles[2], $x_coord[2], $y_coord[2]);
	}

	printf ("%3d %3d %3d ",&product_vs(@hkl1,-1));
	printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[0], $angles[0]+180, -$x_coord[0], -$y_coord[0]);
	printf ("%3d %3d %3d ",&product_vs(@hkl2,-1));
	printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[1], $angles[1]+180, -$x_coord[1], -$y_coord[1]);
	if ($cosines[1] != 0){
		printf ("%3d %3d %3d ",&product_vs(@hkl3,-1));
		printf ("% 4.2f %3d % 3.1f % 3.1f\n", $distance[2], $angles[2]+180, -$x_coord[2], -$y_coord[2]);
	}
	
#Extinction rule
	print ("Shometsusoku (Extinction Rule)\n");
	my @abcreal=&get_abc(@latvec);

#Data: http://www.xtal.iqfr.csic.es/Cristalografia/parte_07_2_1-en.html
#International Tables for X-ray Crystallography kara kiteiru
	print ("Syutsugen jyouken, moshi areba\n");
#@cell=Cubic, Tetragonal, Orthorhombic, Hexagonal/Rhombohedral, Monoclinic
#Cubic/Tetra/Ortho
	if (($abcreal[3] == 90) && ($abcreal[4] == 90) && ($abcreal[5] == 90)){
		if (abs($abcreal[0]-$abcreal[1])>1E-14){
#Ortho
			print ("Orthorhombic:\n");
			print ("Face-centering: h,k,l all odd or all even \n");
			print ("Body-centering: h+k+l=2n \n");
			print ("A-centering: k+l=2n \n");
			print ("C-centering: h+k=2n \n");
			print ("b-glide (100): k=2n @ 0kl\n");
			print ("c-glide (100): l=2n @ 0kl\n");
			print ("n-glide (100): k+l=2n @ 0kl\n");
			print ("d-glide (100): k+l=4n, k=2n, l=2n @ 0kl\n");
			print ("c-glide (010): l=2n @ h0l\n");
			print ("a-glide (010): h=2n @ h0l\n");
			print ("n-glide (010): l+h=2n @ h0l\n");
			print ("d-glide (010): l+h=4n, l=2n, h=2n @ h0l\n");
			print ("a-glide (001): h=2n @ hk0\n");
			print ("b-glide (001): k=2n @ hk0\n");
			print ("n-glide (001): h+k=2n @ hk0\n");
			print ("d-glide (001): h+k=4n, h=2n, k=2n @ hk0\n");
			print ("2_1 screw [100]: h=2n @ h00\n");
			print ("2_1 screw [010]: k=2n @ 0k0\n");
			print ("2_1 screw [001]: l=2n @ 00l\n");
		} else {
			if (abs($abcreal[1]-$abcreal[2])>1E-14){
#Tetra
				print ("Tetragonal:\n");
				print ("Body-centering: h+k+l=2n \n");
				print ("b-glide (100): k=2n @ 0kl\n");
				print ("c-glide (100): l=2n @ 0kl\n");
				print ("n-glide (100): k+l=2n @ 0kl\n");
				print ("d-glide (100): k+l=4n, k=2n, l=2n @ 0kl\n");
				print ("n- and c-glide (1-10): l=2n @ hhl,h-hl\n");
				print ("d-glide (1-10): 2h+l=4n @ hhl,h-hl\n");
				print ("2_1 screw [100]: h=2n @ h00\n");
				print ("4_2 screw [001]: l=2n @ 00l\n");
				print ("4_1 or 4_3 screw [001]: l=4n @ 00l\n");
			} else {
#Cubic
				print ("Cubic:\n");
				print ("Face-centering: h,k,l all odd or all even \n");
				print ("Body-centering: h+k+l=2n \n");
				print ("b-glide (100): k=2n @ 0kl\n");
				print ("c-glide (100): l=2n @ 0kl\n");
				print ("n-glide (100): k+l=2n @ 0kl\n");
				print ("d-glide (100): k+l=4n, k=2n, l=2n @ 0kl\n");
				print ("n- and c-glide (1-10): l=2n @ hhl,h-hl\n");
				print ("d-glide (1-10): 2h+l=4n @ hhl,h-hl\n");
				print ("2_1 or 4_2 screw [100]: h=2n @ h00 \n");
				print ("4_1 or 4_3 screw [100]: h=4n @ h00 \n");
			}
		}
	} elsif (($abcreal[3] == 90) && ($abcreal[4] == 90) && ($abcreal[5] == 120)){
#Hex/Rhombo
		print ("Hexagonal:\n");
		print ("c-glide (11-20): l=2n @ h-h0l\n");
		print ("c-glide (1-100): l=2n @ hh[-2h]l\n");
		print ("6_3 screw [0001]: l=2n @ 000l \n");
		print ("3_1, 3_2, 6_2 or 6_4 screw [0001]: l=3n @ 000l \n");
		print ("6_1 or 6_5 screw [0001]: l=6n @ 000l \n");
		print ("Rhombohedral:\n");
		print ("-h+k+l=3n\n");
		print ("n- and c-glide (1-10): 2h+l=4n @ hhl but no h-hl\n");
	} elsif (($abcreal[3] == 90) && ($abcreal[5] == 90)){
#Monoclinic b-unique C-centering
		print ("Monoclinic:\n");
		print ("C-centering: h+k=2n \n");
		print ("a-glide (001): h=2n @ hk0\n");
		print ("b-glide (001): k=2n @ hk0\n");
		print ("n-glide (001): h+k=2n @ hk0\n");
		print ("2_1 screw [001]: l=2n @ 00l\n");
	}
	exit;
}

sub gcm_direction{
	my @a;
	for (my $i=0;$i<=$#_;$i++){
		push (@a,abs($_[$i])) if ($_[$i] != 0);
	}
	return @_ if ($a[0] eq "");
	my $b=&gcmmany(@a);
	my @c=&product_vs(@_,1/$b);
	for (my $i=0;$i<=$#_;$i++){
		$c[$i]=&round($c[$i]);
	}
	return @c;
}


#Similarity index
#
if ($ARGV[0] eq "--similarity"){
#Henkan gyouretu Determinant, Normalized 2nd invariant, Normalized trace,
# Drift (cartesian x y z), root mean square displacement, normalized root mean square displacement, max displacement, 

# M = STDIN cell [NEW] * (ARGV1 cell [ORIG])^-1

#Subete no gensi zahyou sa ga 0.5 ika ni suru
	&c2d if ($dc eq "C");
	my @arg1cell=&similar($arg2[0]);
#Henkan gyouretu
#M(ARGV1)=(STDIN)
	my @t=&product_mm(@latvec,&invmatrix3(@arg1cell[0..8]));
#determinant
	my $det=&det3(@t)*1;
#2nd
	my $second=$t[0]*$t[0]+$t[4]*$t[4]+$t[8]*$t[8]+2*$t[1]*$t[3]+$t[2]*$t[6]+$t[5]*$t[7];
	$second/=($det**(2/3));
	$second-=3;
#trace
	my $trace=($t[0]+$t[4]+$t[8])/($det**(1/3));
	$trace-=3;
	$det--;
	print ("$det $second $trace ");
#genshi ichi
	for (my $i=0; $i<=$#x; $i++){
		$x[$i]-=$arg1cell[$i*3+9];
		$y[$i]-=$arg1cell[$i*3+10];
		$z[$i]-=$arg1cell[$i*3+11];
	}
#cartesian ni suru
	&d2c;
	my $xsum=&part_sum_array($#x,@x);
	my $ysum=&part_sum_array($#y,@y);
	my $zsum=&part_sum_array($#z,@z);
	my @drift=&product_vs($xsum,$ysum, $zsum, 1/($#x+1));
	print ("@drift ");
	my $sum;
	my $max=0;
	for (my $i=0; $i<=$#x; $i++){
		my $j=sqrt($x[$i]*$x[$i]+$y[$i]*$y[$i]+$z[$i]*$z[$i]);
		$max=$j if ($j>$max);
		$sum+=$j;
	}
	$sum/=($#x+1);
	print ("$sum ");
	$sum/=(((&det3(@latvec))/($#x+1))**(1/3));
	print ("$sum $max \n");
	exit;
}


#Supercell matrix
#
if ($ARGV[0] eq "--supercell_matrix"){
	my @orig_latvec=&get_additional_coordinates($arg2[0],"no_atom_check");
	my @t=&product_mm(@latvec,&invmatrix3(@orig_latvec[0..8]));
	print ("@t \n");
	exit;
}

if ($ARGV[0] eq "--supercell_matrix_round"){
	my @orig_latvec=&get_additional_coordinates($arg2[0],"no_atom_check");
	my @t=&product_mm(@latvec,&invmatrix3(@orig_latvec[0..8]));
	my @t=&round_array(@t);
	print ("@t \n");
	exit;
}


# C) Change POSCAR
# C1) Lattice parameter

if ($ARGV[0] eq "--tidy_basis_vector"){
	my @abcreal=&get_abc(@latvec);
	my @a=&clat(@abcreal,"s");
	@latvec=($a[0],0,0,@a[1..2],0,@a[3..5]);
	&writePOSCAR;
}

if ($ARGV[0] eq "--scale"){
#Direct ni suru
	&c2d if ($dc eq "C");
	if ($arg2[1] eq ""){
		@latvec=&product_mm($arg2[0],0,0,0,$arg2[0],0,0,0,$arg2[0],@latvec);
	} else {
		@latvec=&product_mm($arg2[0],0,0,0,$arg2[1],0,0,0,$arg2[2],@latvec);
	}
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

if ($ARGV[0] eq "--scalexyz"){
#Direct ni suru
	&c2d if ($dc eq "C");
	@latvec=&product_mm(@latvec,$arg2[0],0,0,0,$arg2[1],0,0,0,$arg2[2]);
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

if (($ARGV[0] eq "--strainlattice")||($ARGV[0] eq "--strain_lattice")){
#Direct ni suru
	&c2d if ($dc eq "C");
	$arg2[0]++;
	$arg2[4]++;
	$arg2[8]++;
	@latvec=&product_mm(@latvec,@arg2[0..8]);
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

if ($ARGV[0] eq "--rotate_axis"){
#jiku wo kaiten
	&rotate_axis($arg2[0]);
	&writePOSCAR;
}

if (($ARGV[0] eq "--hklcell_supercell") || ($ARGV[0] eq "--hkl_supercell")){
#jiku wo kaiten
	&get_hkl_cell("supercell",@arg2[0..3]);
	&writePOSCAR;
}

if (($ARGV[0] eq "--hklcell_primitive") || ($ARGV[0] eq "--hkl_primitive")){
#jiku wo kaiten
#c jiku kaeru
	&get_hkl_cell("primitive",@arg2[0..3]);
	&writePOSCAR;
}

#ataerareta hyoumen no nonpolar slab no POSCAR wo syuturyoku 
if (($ARGV[0] eq "--slab_poscar") || ($ARGV[0] eq "--make_slab_poscar")){
	die ("slab atsusa (hikisuu 1) ga ataerarete inai\n") if ($arg2[0] eq "");
	die ("vacuum atsusa (hikisuu 2) ga ataerarete inai\n") if ($arg2[1] eq "");
	die ("H (hikisuu 3) ga ataerarete inai\n") if ($arg2[2] eq "");
	die ("K (hikisuu 4) ga ataerarete inai\n") if ($arg2[3] eq "");
	die ("L (hikisuu 5) ga ataerarete inai\n") if ($arg2[4] eq "");
	&make_surface_slab("normal",@arg2[0..5]);
	exit;
}


#ataerareta hyoumen no nonpolar slab no POSCAR wo H_K_L directory ni syuturyoku
if (($ARGV[0] eq "--slab_poscar_directory") || ($ARGV[0] eq "--make_slab_poscar_directory")){
	&make_surface_slab("directory",@arg2[0..5]);
	die ("slab atsusa (hikisuu 1) ga ataerarete inai\n") if ($arg2[0] eq "");
	die ("vacuum atsusa (hikisuu 2) ga ataerarete inai\n") if ($arg2[1] eq "");
	die ("H (hikisuu 3) ga ataerarete inai\n") if ($arg2[2] eq "");
	die ("K (hikisuu 4) ga ataerarete inai\n") if ($arg2[3] eq "");
	die ("L (hikisuu 5) ga ataerarete inai\n") if ($arg2[4] eq "");
	my $dir=$arg2[2]."_".$arg2[3]."_".$arg2[4];
	print ("WARNING: $dir ga sonzai sinai \n") if (! -e $dir);
	exit;
}


#ataerareta file nohyoumen no nonpolar slab no 
#POSCAR wo H_K_L directory ni syuturyoku
if (($ARGV[0] eq "--slab_poscar_list") || ($ARGV[0] eq "--make_slab_poscar_list")){
	die ("file $arg2[0] ga sonzai sinai\n") if (! -e $arg2[0]);
	die ("slab atsusa (hikisuu 2) ga ataerarete inai\n") if ($arg2[1] eq "");
	die ("vacuum atsusa (hikisuu 3) ga ataerarete inai\n") if ($arg2[2] eq "");
	open LIST, "$arg2[0]";
my $alternative;
	while (my $t=<LIST>){
		my @a=&splitall($t);
		my @b=&make_surface_slab("directory",@arg2[1..2],@a[0..2]);
	}
	close LIST;
	exit;
}


if (($ARGV[0] eq "--slab_step") || ($ARGV[0] eq "--slab_step_poscar") || ($ARGV[0] eq "--slab_facet") || ($ARGV[0] eq "--slab_facet_poscar")) {
	my @a=&make_edge_slab(@arg2);
	exit;
}


if (($ARGV[0] eq "--slab_for_step") || ($ARGV[0] eq "--slab_for_step_poscar ") || ($ARGV[0] eq "--slab_no_facet") || ($ARGV[0] eq "--slab_no_facet_poscar ")){
	$arg2[8]="X" if ($arg2[8] eq "");
	my @a=&make_edge_slab(@arg2[0..2],"XX","XX","XX",@arg2[3..8],"nofacet");
	exit;
}



# C2) Namae kaeru

if (($ARGV[0] eq "--species") || ($ARGV[0] eq "--names")){
	my $count=0;
	for (my $i=0; $i<=$#numspecies; $i++){
		my $name=shift (@arg2);
		for (my $j=0; $j<$numspecies[$i]; $j++){
			$w[$count]=$name;
			$count++;
		}
	}
	&writePOSCAR;
}


if (($ARGV[0] eq "--addnumber") || ($ARGV[0] eq "--add_number")){
	my $count=0;
	for (my $i=0; $i<=$#numspecies; $i++){
		for (my $j=0; $j<$numspecies[$i]; $j++){
			$count++;
			$w[$count-1]=$namespecies[$i].$count;
		}
	}
	&writePOSCAR;
}


if (($ARGV[0] eq "--species6") || ($ARGV[0] eq "--names6")){
	@namespecies=@arg2;
	&writePOSCAR;
}


if (($ARGV[0] eq "--noelement")||($ARGV[0] eq "--no_element")){
	for (my $i=0; $i <= $#namespecies; $i++){
		my @a=split (/_/, $namespecies[$i]);
		$namespecies[$i]=$a[0];
	}
	for (my $i=0; $i <= $#w; $i++){
		$w[$i]=" ";
	}
	&writePOSCAR;
}


if ($ARGV[0] eq "--substitute_atom"){
#substitute atom

	my $xx=$x[$arg2[0]-1];
	my $yy=$y[$arg2[0]-1];
	my $zz=$z[$arg2[0]-1];
	&remove_atoms($arg2[0]);

	push @x,$xx;
	push @y,$yy;
	push @z,$zz;
	push @w,$arg2[1];
	push @namespecies,$arg2[1];
	push @numspecies,1;
	$num_atoms++;
	&writePOSCAR;
	exit;
}


# C3) Gensi zahyo dake

if (($ARGV[0] eq "--incell")||($ARGV[0] eq "--in_cell")){
	&c2d if ($dc eq "C");
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

if (($ARGV[0] eq "--frac") || ($ARGV[0] eq "-f") || ($ARGV[0] eq "-d") || ($ARGV[0] eq "--fract") || ($ARGV[0] eq "--fractional") || ($ARGV[0] eq "--direct")){
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
	&writePOSCAR;
}

if (($ARGV[0] eq "--cart") || ($ARGV[0] eq "-c")){
	if ($dc eq "D"){
		&d2c;
		$dc="C";
	}
	&writePOSCAR;
}


if (($ARGV[0] eq "--selective_true") ||($ARGV[0] eq "--selective_t") ||($ARGV[0] eq "--selective_T")||($ARGV[0] eq "--s_true")||($ARGV[0] eq "--s_t")||($ARGV[0] eq "--s_T")){
	&writePOSCAR_selective("T");
}

if (($ARGV[0] eq "--selective_false") ||($ARGV[0] eq "--selective_f") ||($ARGV[0] eq "--selective_F")||($ARGV[0] eq "--s_false")||($ARGV[0] eq "--s_f")||($ARGV[0] eq "--s_F")){
	&writePOSCAR_selective("F");
}


if ($ARGV[0] eq "--shift"){
#zahyou wo zenbu shift
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&sum_vv($x[$i],$y[$i],$z[$i], @arg2[0..2]);
	}
	&writePOSCAR;
}


if (($ARGV[0] eq "--random") || ($ARGV[0] eq "--randomize")){
#gensi wo zenbu zurasu
	$a=$arg2[0];
	$a=0.001 if ($a eq "");

	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&sum_vv($x[$i],$y[$i],$z[$i], rand(2*$a)-$a, rand(2*$a)-$a, rand(2*$a)-$a);
	}
	&writePOSCAR;
}


if (($ARGV[0] eq "--swap_atom") || ($ARGV[0] eq "--swap_atoms")){
	while ($arg2[1] ne ""){
	#hikisuu check
		die ("gensi #1 ($arg2[0]) ga sonzai sinai\n") if (($x[$arg2[0]-1] eq "") || ($arg2[0] == 0)) ;
		die ("gensi #2 ($arg2[1]) ga sonzai sinai\n") if (($x[$arg2[1]-1] eq "") || ($arg2[1] == 0)) ;
	#swap
		if ($arg2[0] != $arg2[1]){
			($x[$arg2[0]-1],$x[$arg2[1]-1])=($x[$arg2[1]-1], $x[$arg2[0]-1]);
			($y[$arg2[0]-1],$y[$arg2[1]-1])=($y[$arg2[1]-1], $y[$arg2[0]-1]);
			($z[$arg2[0]-1],$z[$arg2[1]-1])=($z[$arg2[1]-1], $z[$arg2[0]-1]);
			($w[$arg2[0]-1],$w[$arg2[1]-1])=($w[$arg2[1]-1], $w[$arg2[0]-1]);
		}
		shift @arg2;
		shift @arg2;
	}
	&writePOSCAR;
}



if ($ARGV[0] eq "--swap_species"){
#hikisuu check
	die ("genso #1 ($arg2[0]) ga sonzai sinai\n") if (($numspecies[$arg2[0]-1] eq "") || ($arg2[0] == 0)) ;
	die ("genso #2 ($arg2[1]) ga sonzai sinai\n") if (($numspecies[$arg2[1]-1] eq "") || ($arg2[1] == 0)) ;
#$arg2[0]<$arg2[1] ni suru 
	@arg2[0..1]=($arg2[1], $arg2[0]) if ($arg2[0]>$arg2[1]);
#genso no jyunban wo kaeru
	my @orig_order=(0..$#numspecies);
	my @new_order=@orig_order;
	$new_order[$arg2[0]-1]=$orig_order[$arg2[1]-1];
	$new_order[$arg2[1]-1]=$orig_order[$arg2[0]-1];
#numspecies
	my @numspecies_new=@numspecies[@new_order];
#namespecies
	@namespecies=@namespecies[@new_order];
#swap a lot of sites, note $arg2[0]<$arg2[1]
	@x=&swap_species_coordinates ($arg2[0]-1,$arg2[1]-1,@numspecies_new,@x);
	@y=&swap_species_coordinates ($arg2[0]-1,$arg2[1]-1,@numspecies_new,@y);
	@z=&swap_species_coordinates ($arg2[0]-1,$arg2[1]-1,@numspecies_new,@z);
	@w=&swap_species_coordinates ($arg2[0]-1,$arg2[1]-1,@numspecies_new,@w);
	@numspecies=@numspecies_new;
	&writePOSCAR;
}

sub swap_species_coordinates{
# genso $_[0] and $_[1] ($_[0] < $_[1]) no zahyo wo koukan
# atarasii genso no jyunban ha @_[2..($#numspecies+1)];
# zahyo wa @_[($#numspecies+2)..]
	my @a=@_;
#genso data
	my $site1=shift @a;
	my $site2=shift @a;
	my @numspecies_new=splice(@a,0,$#numspecies+1);

	my @tmp1;
	my @tmp2;
#soto kara kesu
	@tmp1=splice(@a,&part_sum_array($site2,@numspecies),$numspecies[$site2]);
	@tmp2=splice(@a,&part_sum_array($site1,@numspecies),$numspecies[$site1]);
#naka kara tasu
	splice(@a,&part_sum_array($site1,@numspecies_new),0,@tmp1);
	splice(@a,&part_sum_array($site2,@numspecies_new),0,@tmp2);
#return
	return (@a);
}

if ($ARGV[0] eq "--similar"){
#Subete no gensi zahyou sa ga 0.5 ika ni suru
	&c2d if ($dc eq "C");
	&similar($ARGV[1]);
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

#similar: POSCAR wo $_[0] ni chikadukeru
#@ref wo kaesu

sub similar{
	my @ref=&get_additional_coordinates($_[0]);
	for (my $i=0;$i<$num_atoms;$i++){
		while (($x[$i]-$ref[$i*3+9]) <= -0.5){
			$x[$i]++;
		}
		while (($x[$i]-$ref[$i*3+9]) > 0.5){
			$x[$i]--;
		}
		while (($y[$i]-$ref[$i*3+10]) <= -0.5){
			$y[$i]++;
		}
		while (($y[$i]-$ref[$i*3+10]) > 0.5){
			$y[$i]--;
		}
		while (($z[$i]-$ref[$i*3+11]) <= -0.5){
			$z[$i]++;
		}
		while (($z[$i]-$ref[$i*3+11]) > 0.5){
			$z[$i]--;
		}
	}
	return (@ref);
}



if (($ARGV[0] eq "--nebimages")||($ARGV[0] eq "--neb_images")){
#$arg2[1]=totyu image no kazu
	die ("neb image no kazu ga okasii: $arg2[1] \n") if ($arg2[1] !~ /^\d{1,2}$/);	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
	my @ref=&similar($arg2[0]);
#@x,@y,@z: saisyo copy to @startx,@starty,@startz
#@ref: saigo copy to @endx,@endy,@endz
	my @startx=@x;
	my @starty=@y;
	my @startz=@z;
	my (@endx, @endy, @endz);
	splice (@ref, 0, 9);
	for (my $i=0; $i<$num_atoms; $i++){
		$endx[$i]=$ref[$i*3];
		$endy[$i]=$ref[$i*3+1];
		$endz[$i]=$ref[$i*3+2];
	}
	my $maxdir=$arg2[1]+1;
	for (my $j=0;$j<=$maxdir;$j++){
		my $dir;
		if ($j<10){
			$dir="0".$j;
		} else {
			$dir=$j;
		}
		for (my $i=0; $i<$num_atoms; $i++){
			$x[$i]=$startx[$i]*($j/$maxdir)+$endx[$i]*(1-$j/$maxdir);
			$y[$i]=$starty[$i]*($j/$maxdir)+$endy[$i]*(1-$j/$maxdir);
			$z[$i]=$startz[$i]*($j/$maxdir)+$endz[$i]*(1-$j/$maxdir);
		}
		system ("mkdir $dir");
		&writePOSCAR ("$dir/POSCAR");
	}
	exit;
}

#symmetrize
if ($ARGV[0] eq "--symmetrize"){
	die ("hikisuu ga 12 nai \n") if ($arg2[11] eq "");
#cutoff
	my $dist=0.3;
	if ($arg2[12] ne ""){
		$dist=$arg2[12]*1 ;
		die ("dist $arg2[12] ga kazu de nai?") if ($dist <= 0);
	}
	&symmetrize (@arg2[0..11],$dist);
	&writePOSCAR;
	exit
}


# C4) Kousi henkei

if ($ARGV[0] eq "--niggli"){
#Direct ni suru
	&c2d if ($dc eq "C");
	my @latvec_niggli=&get_niggli_latvec(@latvec);
	my @conversion_matrix=&round_array(&product_mm(@latvec,&invmatrix3(@latvec_niggli)));
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i], @conversion_matrix);
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	@latvec=@latvec_niggli;
	&d2c if ($dc eq "C");
	&writePOSCAR;
}


if (($ARGV[0] eq "--reciniggli")||($ARGV[0] eq "--reci_niggli")){
#Direct ni suru
	&c2d if ($dc eq "C");
#Reci
	my @latvec_reci=&tmatrix3(&invmatrix3(@latvec));
	my @latvec_reci_niggli=&get_niggli_latvec(@latvec_reci);
	my @latvec_niggli=&tmatrix3(&invmatrix3(@latvec_reci_niggli));
	my @conversion_matrix=&round_array(&product_mm(@latvec,&invmatrix3(@latvec_niggli)));
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i], @conversion_matrix);		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	@latvec=@latvec_niggli;
	&d2c if ($dc eq "C");
	&writePOSCAR;
}


if (($ARGV[0] eq "--maxortho")||($ARGV[0] eq "--max_ortho")){
#Algorithm Based on "A 3-Dimensional Lattice Reduction Algorithm"
#(Igor Semaev, Cryptography and Lattices. Volume 2146 of the series 
#Lecture Notes in Computer Science pp 181-193) p188-p189 
	&d2c if ($dc eq "D");
	&maxortho;
	&c2d;
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

#shinkuusou atsusa kaeru
if (($ARGV[0] eq "--change_slab_vacuum_thickness")||($ARGV[0] eq "--change_slab_vacuum_thickness_ortho")){
	die ("atsusa ga hikisuu ni nai \n") if ($arg2[0] eq "");
	&c2d if ($dc eq "C");
#ichiou site_in_cell
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
#center wo z=0
	my ($minz, $maxz)=(&min(@z), &max(@z));
	my $centerz=($minz+$maxz)*0.5;
	for (my $i=0; $i<$num_atoms; $i++){
		$z[$i]-=$centerz;
	}
#c wo kaeru
	my @ab=&cross_vv(@latvec[0..5]);
	my $vol=&product_vv(@ab,@latvec[6..8]);
	my $abnorm=&norm(@ab);
	my $height=$vol/$abnorm;
	my $slab=($maxz-$minz)*$height;
	my @reduce=&product_vs(@ab,($height-$slab-$arg2[0])/$abnorm);
	&d2c;
	if ($ARGV[0] eq "--change_slab_vacuum_thickness"){
#legacy definition
		@latvec[6..8]=&diff_vv(@latvec[6..8],@reduce);
	} elsif ($ARGV[0] eq "--change_slab_vacuum_thickness_ortho"){
		@latvec[6..8]=&product_vs(@ab,($slab+$arg2[0])/$abnorm);
	} else {
		die ("change_slab_vacuum_thickness: arienai\n");
	}
	&c2d;
	for (my $i=0; $i<$num_atoms; $i++){
		$z[$i]+=0.5;
	}
	&d2c if ($dc eq "C");
	&writePOSCAR;
}




# C5) Gensi kezuru

#unique sites
if ($ARGV[0] eq "--unique"){
	&c2d if ($dc eq "C");
	#in_cell first
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	&unique($arg2[0]);
	&d2c if ($dc eq "C");
	&writePOSCAR;
}

#slab wo tsukuru
if (($ARGV[0] eq "--makeslab")||($ARGV[0] eq "--make_slab")){
	&make_slab(@arg2[0..3]);
	&writePOSCAR;
}


if (($ARGV[0] eq "--rmatom")||($ARGV[0] eq "--rm_atom")||($ARGV[0] eq "--rmatoms")||($ARGV[0] eq "--rm_atoms")){
	my @rmlist=&expand_array(@arg2);
	&remove_atoms(@rmlist);
	&writePOSCAR;
}

if (($ARGV[0] eq "--rmspecies")||($ARGV[0] eq "--rm_species")){
	my @rmsitelist=&expand_array(@arg2);
	my @rmatomlist;
	for (my $i=1; $i<=$num_atoms; $i++){
		my @a=&nameatom($i);
		for (my $j=0; $j<=$#rmsitelist;$j++){
			push (@rmatomlist, $i) if ($a[1] == $rmsitelist[$j]);
		}
	}
	&remove_atoms(@rmatomlist);
	if ($x[0] eq ""){
		print ("--rm_species no kekka gensi ga sonzai sinai: no POSCAR \n");
		exit;
	}
	&writePOSCAR;
}



#daruma otoshi no youni gensi wo kezuru
if ($ARGV[0] eq "--daruma") {
#gensi 2 (x2,y2,z2) wo genshi 1 (x1,y1,z1) no ichi ni mottekuru 
#z>=z2 no gensi wo (x1-x2,y1-y2,z1-z2) idou
#z2>z>=z1 no gensi ha sakujyo 

	&c2d if ($dc eq "C");
	my @t=($x[$arg2[0]-1]-$x[$arg2[1]-1],$y[$arg2[0]-1]-$y[$arg2[1]-1],$z[$arg2[0]-1]-$z[$arg2[1]-1]);
	my $zmax=$z[$arg2[1]-1];
	my $zmin=$z[$arg2[0]-1];
	my @rmlist;
	for (my $i=0; $i<$num_atoms; $i++){
		if ($z[$i]>=$zmax){
			($x[$i],$y[$i],$z[$i])=&sum_vv($x[$i],$y[$i],$z[$i],@t);
		} elsif ($z[$i] >= $zmin){
			push (@rmlist,$i+1);
		}
	}
	&remove_atoms(@rmlist);
	&d2c if ($dc eq "C");
	&writePOSCAR;
}


sub remove_atoms{
#hikisuu no gensi wo kesu
#hikisuuno saisyo no gensi no bango ha 1 tosuru
	my @rmlist = sort { $a <=> $b } @_;
#ookii kazu kara kesu
	for (my $i=$#rmlist; $i>=0;$i--){
		splice(@x,$rmlist[$i]-1,1);
		splice(@y,$rmlist[$i]-1,1);
		splice(@z,$rmlist[$i]-1,1);
		splice(@w,$rmlist[$i]-1,1);
		my @a=&nameatom($rmlist[$i]);
		$numspecies[$a[1]-1]--;
		$num_atoms--;
	}
	&remove_zero_atom_species;
}


if (($ARGV[0] eq "--supercell")||($ARGV[0] eq "--transform")||($ARGV[0] eq "--change_of_basis")||($ARGV[0] eq "--cob")){
	&get_supercell(@arg2);
	&writePOSCAR;
}


if (($ARGV[0] eq "--addatom")||($ARGV[0] eq "--add_atom")){
	$arg2[3]="X" if ($arg2[3] eq "");
	push @x,$arg2[0];
	push @y,$arg2[1];
	push @z,$arg2[2];
	push @w,$arg2[4];
	push @namespecies,$arg2[3];
	push @numspecies,1;
	$num_atoms++;
	&writePOSCAR;
}


if (($ARGV[0] eq "--addatomexist")||($ARGV[0] eq "--add_atom_exist")){
	#species ga aru ka?
	die ("species bangou $arg2[3] ga seisuu de nai \n") if ($arg2[3] !~ /^[0-9]+$/);
	my $max_species=$#numspecies+1;
	die ("species bangou $arg2[3] ga sonzai sinai; 1 to $max_species kara erabu") if (! &within($arg2[3],1,$max_species));
	splice(@x,&part_sum_array($arg2[3]-1,@numspecies),0,$arg2[0]);
	splice(@y,&part_sum_array($arg2[3]-1,@numspecies),0,$arg2[1]);
	splice(@z,&part_sum_array($arg2[3]-1,@numspecies),0,$arg2[2]);
	splice(@w,&part_sum_array($arg2[3]-1,@numspecies),0,$arg2[4]);
	$numspecies[$arg2[3]-1]++;
	$num_atoms++;
	&writePOSCAR;
}

if (($ARGV[0] eq "--superimpose")||($ARGV[0] eq "--superimpose_poscar")){
#superimpose poscar
	my $count=1;
	my $current_file=shift(@arg2);
	while ($current_file ne ""){
		die ("superimpose: no file $current_file \n") if (! -e $current_file);
		$count++;
		my @ref=&get_additional_coordinates($current_file);
		splice (@ref, 0, 9);
		for (my $i=0; $i<$num_atoms; $i++){
			splice(@x,($i+1)*$count-1,0,$ref[$i*3]);
			splice(@y,($i+1)*$count-1,0,$ref[$i*3+1]);
			splice(@z,($i+1)*$count-1,0,$ref[$i*3+2]);
		}
		$current_file=shift(@arg2);
	}
	$num_atoms*=$count;
	@numspecies=&product_vs(@numspecies,$count);
	&writePOSCAR;
}

#recover slab: slab no nokori wo isometry de fukugen
if (($ARGV[0] eq "--recover_slab")||($ARGV[0] eq "--recoverslab")){
	my (@x2, @y2, @z2, @w2);
	my (@x3, @y3, @z3, @w3);
	for (my $i=0; $i<$num_atoms;$i++){
		my @rotate=&product_mv(@arg2[0..8],$x[$i],$y[$i],$z[$i]);
		($x2[$i],$y2[$i],$z2[$i])=&sum_vv(@rotate,@arg2[9..11]);
		$w2[$i]=$w[$i];
	}
	for (my $i=0; $i<=$#numspecies;$i++){
		push (@x3,splice(@x2,0,$numspecies[$i]));
		push (@x3,splice(@x,0,$numspecies[$i]));
		push (@y3,splice(@y2,0,$numspecies[$i]));
		push (@y3,splice(@y,0,$numspecies[$i]));
		push (@z3,splice(@z2,0,$numspecies[$i]));
		push (@z3,splice(@z,0,$numspecies[$i]));
		push (@w3,splice(@w2,0,$numspecies[$i]));
		push (@w3,splice(@w,0,$numspecies[$i]));
	}
	@numspecies=&product_vs(@numspecies,2);
	$num_atoms*=2;
	for (my $i=0; $i<$num_atoms;$i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x3[$i],$y3[$i],$z3[$i]);
	}
	&unique($arg2[12]);
	&writePOSCAR;
}


#USPEX slab wo tsukuru
if (($ARGV[0] eq "--uspex_slab")||($ARGV[0] eq "--uspexslab")){
	my @numspecies_orig=@numspecies;
#inversion @ genten or mirror @ z=0
	&get_primitive_for_surface_reconstruction(@arg2[0..2],$arg2[4]);
	my @filename=("POSCAR.INV.A.1.","POSCAR.INV.A.2.","POSCAR.INV.B.1.","POSCAR.INV.B.2.","POSCAR.MIR.A.1.","POSCAR.MIR.A.2.","POSCAR.MIR.B.1.","POSCAR.MIR.B.2.");
#supercell, USPEX cell wo tsukuru
	foreach my $i (@filename){
		my $j=$i."primitive";
		if (-e $j){
		print ("$j wo syori\n");
#poscar wo yomu
			my @poscar_info=&get_additional_coordinates($j);
			my $line7=`head -7 $j | tail -n 1`;
			@numspecies=&splitall($line7);
			$num_atoms=&numatoms;
			@latvec=splice @poscar_info, 0, 9;
			@x=@y=@z=@w=();
			for (my $k=0; $k < $num_atoms; $k++){
				$x[$k]=shift @poscar_info;
				$y[$k]=shift @poscar_info;
				$z[$k]=shift @poscar_info;
				$w[$k]=" ";
			}
#supercell size wo kimeru
			my @ab=&cross_vv(@latvec[0..5]);
			my $vol=&product_vv(@ab,@latvec[6..8])*$scale*$scale*$scale;
			my $abnorm=&norm(@ab)*$scale*$scale;
			my $height=$vol/$abnorm;
			my $n=&ceiling($arg2[3]/$height);
			&get_supercell(1,1,$n);
			&d2c;
			&maxortho;
			&c2d;
#supercell wo kaku
			my $outfile=$i."supercell";
		
			&writePOSCAR($outfile);
#koko kara USPEX version!
			my @latpar=&get_abc(@latvec);
			my $newheight=$height*$n;
			my $maxz=&max(@z);
			my @lat=&clat(@latpar, "s");
			@latvec=($lat[0],0,0,@lat[1..2],0,@lat[3..5]);
			$outfile=$i."supercell_rotated";
			&writePOSCAR($outfile);

			&d2c;
			@latvec=($lat[0],0,0,@lat[1..2],0,0,0,$newheight*$maxz*1.0001);
			$outfile=$i."USPEX";
			&c2d;
			&writePOSCAR($outfile);
#gensi suu ga kowasareteiru kamo sirenai : naosu
			@numspecies=@numspecies_orig;
		}
	}
	exit;
}


if ($ARGV[0] eq "--addimages"){
#Direct ni suru
	&c2d if ($dc eq "C");
#tyuokuhoutai cell wo tukuru
	my @latvec_orig=@latvec;
	die ("--addimages; minx >= maxx, die \n") if ($arg2[0] >= $arg2[1]);
	die ("--addimages; miny >= maxy, die \n") if ($arg2[2] >= $arg2[3]);
	die ("--addimages; minz >= maxz, die \n") if ($arg2[4] >= $arg2[5]);
	my $minx=&floor($arg2[0]);
	my $miny=&floor($arg2[2]);
	my $minz=&floor($arg2[4]);
	my $maxx=&ceiling($arg2[1]+0.00001);
	my $maxy=&ceiling($arg2[3]+0.00001);
	my $maxz=&ceiling($arg2[5]+0.00001);
	my (@xnew, @ynew, @znew, @wnew);
	my $num_atoms_new;
	my @numspecies_new=&diff_vv(@numspecies,@numspecies);
	
	for (my $i=0;$i<$num_atoms;$i++){
		for (my $j=$minx;$j<$maxx;$j++){
			for (my $k=$miny;$k<$maxy;$k++){
				for (my $l=$minz;$l<$maxz;$l++){
					if (&within($x[$i]+$j,$arg2[0],$arg2[1])){
						if (&within($y[$i]+$k,$arg2[2],$arg2[3])){
							if (&within($z[$i]+$l,$arg2[4],$arg2[5])){
								push (@xnew,$x[$i]+$j);
								push (@ynew,$y[$i]+$k);
								push (@znew,$z[$i]+$l);
								push (@wnew,$w[$i]);
								my @t=&nameatom($i+1);
								$numspecies_new[$t[1]-1]++;
								$num_atoms_new++;
							}
						}
					}
				}
			}
		}
	}
	@x=@xnew;
	@y=@ynew;
	@z=@znew;
	@w=@wnew;
	@numspecies=@numspecies_new;
	$num_atoms=$num_atoms_new;
	&remove_zero_atom_species;
	&d2c if ($dc eq "C");
#gensisuu wo update
	&writePOSCAR;
}

# D) Gaibu interface

if (($ARGV[0] eq "--cif") || ($ARGV[0] eq "--cif6")) {
	&writeCIF;
}

if (($ARGV[0] eq "--xyz") || ($ARGV[0] eq "--xyz6")) {
	&writexyz;
}


#E) Standard conventional/primtive/band zu

#Kessyougaku conventional
if (($ARGV[0] eq "--conventional") || ($ARGV[0] eq "--bposcar") || ($ARGV[0] eq "--crystallographic_conventional")|| ($ARGV[0] eq "--cc")) {
	&get_crystallographic_cell("conventional");
	&writePOSCAR;
}

#Kessyougaku primitive
if (($ARGV[0] eq "--crystallographic_primitive") || ($ARGV[0] eq "--cp")) {
	&get_crystallographic_cell("primitive");
	&writePOSCAR;
}

#Hinuma reduced (triclinic only)
if ($ARGV[0] eq "--hinuma_reduced") {
	&get_hinuma_reduced("primitive");
	&writePOSCAR;
}

#Standard conventional
if (($ARGV[0] eq "--standard_conventional") || ($ARGV[0] eq "--sc")) {
	&get_standard_SC("conventional");
	&writePOSCAR;
}

#Standard primitive
if (($ARGV[0] eq "--standard_primitive") || ($ARGV[0] eq "--sp")) {
	&get_standard_SC("primitive");
	&writePOSCAR;
}


#BZ point kanren
# BZ koutaisyouten list
if (($ARGV[0] eq "--bzpoint") || ($ARGV[0] eq "--bzpoint_notrs")){
	my @points;
	if ($ARGV[0] eq "--bzpoint"){
		$arg2[0]="KPOINTS.Cbzpoint" if ($arg2[0] eq "");
		@points=&bz_point("Y");
	} else {
		$arg2[0]="KPOINTS.Cbzpoint_notrs" if ($arg2[0] eq "");
		@points=&bz_point;
	}
	shift @points;
	open OUT, ">$arg2[0]";
	for (my $i=0; 4*$i < $#points ; $i++){
		printf OUT ("% 2.7f % 2.7f % 2.7f 0 ",@points[$i*4..$i*4+2]);
		print OUT ("! $points[$i*4+3] \n");
	}
	close OUT;
	&writePOSCAR;
}

if ($ARGV[0] eq "--scbzpoint"){
	$arg2[0]="KPOINTS.SCbzpoint" if ($arg2[0] eq "");
	my @points=&bz_point_SC;
	shift @points;
	open OUT, ">$arg2[0]";
	for (my $i=0; 4*$i < $#points ; $i++){
		printf OUT ("% 2.7f % 2.7f % 2.7f 0 ",@points[$i*4..$i*4+2]);
		print OUT ("! $points[$i*4+3] \n");
	}
	close OUT;
	&writePOSCAR;
}


# BZ koutaisyouten kpf
if (($ARGV[0] eq "--bzpointkpf") || ($ARGV[0] eq "--bzpointkpf_notrs")){
	my @points;
	if ($ARGV[0] eq "--bzpointkpf"){
		$arg2[0]="Cbzpoint.kpf" if ($arg2[0] eq "");
		@points=&bz_point("Y");
	} else {
		$arg2[0]="Cbzpoint_notrs.kpf" if ($arg2[0] eq "");
		@points=&bz_point;
	}
	shift @points;
	open OUT, ">$arg2[0]";
	print OUT ("Real form of k-point coordinates (kx,ky,kz,label):\n");
	for (my $i=0; 4*$i < $#points ; $i++){
		printf OUT (" % 2.7f % 2.7f % 2.7f K.%d %s \n",@points[$i*4..$i*4+2],$i+1,$points[$i*4+3]);
	}
	close OUT;
	&writePOSCAR;
}


if ($ARGV[0] eq "--scbzpointkpf"){
	$arg2[0]="SCbzpoint.kpf" if ($arg2[0] eq "");
	my @points=&bz_point_SC;
	shift @points;
	open OUT, ">$arg2[0]";
	print OUT ("Real form of k-point coordinates (kx,ky,kz,label):\n");
	for (my $i=0; 4*$i < $#points ; $i++){
		printf OUT (" % 2.7f % 2.7f % 2.7f K.%d %s \n",@points[$i*4..$i*4+2],$i+1,$points[$i*4+3]);
	}
	close OUT;
	&writePOSCAR;
}


# BZ kpath
if (($ARGV[0] eq "--kpath") || ($ARGV[0] eq "--kpath_notrs")){
	$arg2[0]=0.025 if ($arg2[0] eq "");
	if ($ARGV[0] eq "--kpath"){
		$arg2[1]="KPOINTS.Ckpath_".${arg2[0]} if ($arg2[1] eq "");
		&get_kpath_kpoints($arg2[0],"Cryst", "Y");
	} else {
		$arg2[1]="KPOINTS.Ckpath_notrs_".${arg2[0]} if ($arg2[1] eq "");
		&get_kpath_kpoints($arg2[0],"Cryst");
	}
	system ("mv tsubo_temp_get_kpath $arg2[1]");
	&writePOSCAR;
}




if ($ARGV[0] eq "--sckpath"){
	$arg2[0]=0.025 if ($arg2[0] eq "");
	$arg2[1]="KPOINTS.SCkpath_".${arg2[0]} if ($arg2[1] eq "");
	&get_kpath_kpoints($arg2[0],"SC");
	system ("mv tsubo_temp_get_kpath $arg2[1]");
	&writePOSCAR;
}


# BZ kpath kpf
if (($ARGV[0] eq "--kpathkpf") || ($ARGV[0] eq "--kpathkpf_notrs")){
	if ($ARGV[0] eq "--kpathkpf"){
		$arg2[0]="defaultCkpath.kpf" if ($arg2[0] eq "");
		&get_kpath_kpf("Y");
	} else {
		$arg2[0]="defaultCkpath_notrs.kpf" if ($arg2[0] eq "");
		&get_kpath_kpf;
	}
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}

if ($ARGV[0] eq "--sckpathkpf"){
	$arg2[0]="defaultSCkpath.kpf" if ($arg2[0] eq "");
	&get_kpath_kpf_SC;
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}

# BZ kpath phonopy
if ($ARGV[0] eq "--kpathphonopy") {
	$arg2[0]="defaultCkpath.phonopy" if ($arg2[0] eq "");
	&get_kpath_phonopy("Cryst","Y");
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}

if ($ARGV[0] eq "--sckpathphonopy"){
	$arg2[0]="defaultSCkpath.phonopy" if ($arg2[0] eq "");
	&get_kpath_phonopy("SC");
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}


# BZ kpath band
if (($ARGV[0] eq "--kpath_band") || ($ARGV[0] eq "--kpath_band_notrs")){
	if ($ARGV[0] eq "--kpath_band"){
		$arg2[0]="KPOINTS.Ckpath_band" if ($arg2[0] eq "");
		&get_kpath_kpoints(0.025, "Cryst", "Y");
	} else {
		$arg2[0]="KPOINTS.Ckpath_band_notrs" if ($arg2[0] eq "");
		&get_kpath_kpoints(0.025, "Cryst");
	}
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}

if ($ARGV[0] eq "--sckpath_band"){
	$arg2[0]="KPOINTS.SCkpath_band" if ($arg2[0] eq "");
	&get_kpath_kpoints(0.025, "SC");
	system ("mv tsubo_temp_get_kpath $arg2[0]");
	&writePOSCAR;
}


# BZ kpath mass calculation
if (($ARGV[0] eq "--kpath_mass") || ($ARGV[0] eq "--kpath_mass_notrs")){
	if ($ARGV[0] eq "--kpath_mass") {
		$arg2[0]="KPOINTS.Ckpath_mass" if ($arg2[0] eq "");
		&get_kpath_mass(0.002, "Cryst","Y");
	} else {
		$arg2[0]="KPOINTS.Ckpath_mass_notrs" if ($arg2[0] eq "");
		&get_kpath_mass(0.002, "Cryst");
	}
	system ("mv tsubo_temp1 $arg2[0]");
	&writePOSCAR;
}

if ($ARGV[0] eq "--sckpath_mass"){
	$arg2[0]="KPOINTS.SCkpath_mass" if ($arg2[0] eq "");
	&get_kpath_mass(0.002, "SC");
	system ("mv tsubo_temp1 $arg2[0]");
	&writePOSCAR;
}

	print <<"EOM";
Tadashii option ga arimasen
-h wo tameshite kudasai
Moshika shite..

tsubo.pl --data < POSCAR 
tsubo.pl --cif6 < POSCAR 
tsubo.pl --dist D < POSCAR 
tsubo.pl --sgsymbol < POSCAR 

tsubo.pl --direct < POSCAR 
tsubo.pl --cart < POSCAR 
tsubo.pl --cc < POSCAR 
tsubo.pl --cp < POSCAR 
tsubo.pl --niggli < POSCAR 
tsubo.pl --supercell P Q R S T U V W X < POSCAR 

tsubo.pl --kpoints Spacins < POSCAR 
tsubo.pl --kpoints_slab Spacins < POSCAR 

tsubo.pl --kpath_band_notrs < POSCAR 
tsubo.pl --kpath_mass_notrs < POSCAR 

tsubo.pl --termination_polarity H K L [asis] 
tsubo.pl --hkl_supercell H K L [asis] 
tsubo.pl --hkl_primitive H K L [asis] 
tsubo.pl --slab_poscar slab_thickness vacuum_thickness H K L [asis] < POSCAR 

tsubo.pl --unique_nonpolar H K L cap < POSCAR > tempfile1
tsubo.pl --termination_polarity_list tempfile1 < POSCAR > tempfile2
egrep 'A|B' tempfile2 | head -N > tempfile3
tsubo.pl --slab_poscar_list tempfile3 slab_thickness vacuum_thickness <POSCAR

tsubo.pl --similar POSCAR1 <POSCAR2
tsubo.pl --nebimages POSCAR1 N <POSCAR2

EOM

die ("\n");

#### Main syuuryou



sub help{
	print <<"EOM";
TSUBOGURUMA by Hiroaki Nissho

Kihonteki tsukaikata: tsubo.pl --command <POSCAR

Commands:

-h -help --help : Help wo hyouji
--version : Version wo hyouji

--clat : a,b,c,alpha,beta,gamma kara lattice vector wo tsukuru (POSCAR iranai)
--symop_info --symmetry_operation_info : taisyou sousa no jyouhou wo dasu

--data : kousi teisuu nado wo hyouji
--data_full : kousi teisuu nado wo hyouji (syousuuten 10 keta made)
--dataniggli : niggli reduced cell no kousi teisuu nado wo hyouji
--natoms --numatoms : Gensi no kazu wo hyouji
--nspecies --numspecies : Genso no kazu wo hyouji
--stoichiometry : soseihi wo kaesu
--which_atom_number : chikai gensi no bangou wo kaesu
--dist D : Kyori D ika no gensi wo hyouji
--dist_large D : Kyori D ika no gensi wo hyouji (kensaku hanni wo hirogeru)
--dist_min : Saisyou genshikan kyori wo hyouji
--dist_atom_list D (list) : (list) no gensi ni tsuite kyori D ika no gensi wo hyouji
--solid_angle d : hyoumen site jidou hantei
--coordination [Cutoff1] [Cutoff2] : haiisuu wo hyouji
--coordination_large [Cutoff1] [Cutoff2]  : haiisuu wo hyouji (kensaku hanni wo hirogeru)
--maxgap : c jiku houkou max gap wo hyouji
--volperatom : Volume/atom wo hyouji
--nameatom N : gensi N no gensomei (6 gyome) to genso bango (7 gyome) wo kaesu
--kpoints Spacing [Mesh] : KPOINTS wo tsukuru
--kpoints_even Spacing [Mesh] : even KPOINTS wo tsukuru
--kpoints_slab Spacing [Mesh] : KPOINTS wo tsukuru (slab model)
--kpoints_even_slab Spacing [Mesh] : even KPOINTS wo tsukuru (slab model)
--slabarea : slab area wo kaesu
--average_atom N1 N2.. : gensi no heikin zahyou wo kaesu
--sgnumber : Kuukan gun bangou wo kaesu
--sgsymbol : Kuukan gun symbol wo kaesu
--2Dsgnumber z : 3D poscar no z danmen no 2D Kuukan gun bangou wo kaesu
--2Dsgsymbol z : 3D poscar no z danmen no 2D Kuukan gun symbol wo kaesu
--pgsymbol : Ten gun symbol wo kaesu
--crystalsystem : Kesshou kei wo kaesu (Triclinic nado)
--bravais : Bravais lattice wo kaesu (aP nado)
--sclatticetype : Kousi no shurui wo kaesu (eg CUB, ORCF2)
--normal_vector H K L length : HKL hyoumen ni suichoku de nagasa length no vector wo hyouji
--isometry [info] : Subete no isometry wo hyouji
--isometry_nonpolar [info]: nonpolar isometry wo hyouji
--isometry_nonpolar_square : nonpolar isometry da ga  nijyou site identity  no isometry wo hyouji 
--surface_isometry H K L [asis] : HKL hyoumen no slab no taisyousei wo
hyouji
--unique_nonpolar [h k l] [cap] : nonpolar surface wo chiisai jun ni kaesu
--to_unique_nonpolar H K L : unique nonpolar orientation wo kaesu
--isometry_nonpolar_image_file file w11 w12 w13 w21 w22 w23 w31 w32 w33 w1 w2 w3 : isometry ni yoru file no image wo tsukuru
--termination_polarity H K L [asis] : HKL hyoumen no slab no osusume wo hyouji
--termination_polarity_full H K L [asis] : HKL hyoumen no slab no polarity wo hyouji
--termination_polarity_list file : file no hyoumen no slab no polarity wo hyouji
--inplane_lattice_vector H K L [asis] : HKL hyoumen no inplane lattice vector wo hyouji
--terrace_vicinal_vector H K L He Ke Le [asis] : HKL hyoumen de edge ga He Ke Le houkou no terrace to vicinal vector wo hyouji
--orientation_relative_all H K L Hmax Kmax Lmax : HKL hyoumen ni taisuru hoka no orientation wo hyouji
--relative_orientation_facet H K L Hmax Kmax Lmax : HKL hyoumen de facet ga tsukuresou na orientation wo hyouji 
--facet_surface_energy H K L Hmax Kmax Lmax library : facet surface no energy wo hyouji 
--stepenergy Esurf_vicinal Esurf_terrace N1 N2 ... : step energy wo keisan
--equivalent_site : Touka na site wo hyouji
--2ddiffraction H K L : H K L houkou no 2D diffraction pattern wo kaesu
--similarity POSCAR1 : Similarity index wo kaesu
--supercell_matrix original_POSCAR < new_POSCAR  : supercell matrix wo kaesu
--supercell_matrix_round original_POSCAR < new_POSCAR  : supercell matrix wo kaesu

--tidy_basis_vector : basis vector wo kirei ni suru
--scale A : POSCAR no kousi teisuu wo A bai suru
--scale A B C: POSCAR no kousi teisuu wo sorezore A,B,C bai suru
--scalexyz A B C: POSCAR no x,y,z zahyo wo A,B,C bai suru
--strainlattice A B C D E F G H I:kousi wo hizumaseru
--rotate_axis abc: POSCAR no kitei vector wo mawasu
--hklcell_supercell --hkl_supercell H K L [asis] : HKL houkou no supercell wo tsukuru
--hklcell_primitive --hkl_primitive H K L [asis] : HKL houkou no primitive cell wo tsukuru
--slab_poscar --make_slab_poscar slab_thickness vacuum_thickness H K L [asis] : HKL hyoumen no slab wo kaesu
--slab_poscar_directory --make_slab_poscar_directory slab_thickness vacuum_thickness H K L [asis] : directory wo tsukutte HKL hyoumen no slab wo kaesu
--slab_poscar_list --make_slab_poscar_list slab_thickness vacuum_thickness file : directory wo tsukutte file no hyoumen no slab wo kaesu
--slab_step he ke le ht kt lt hv kv lv slab_thickness vacuum_thickness [asis] : step ga aru model wo tsukuru
--slab_for_step he ke le hv kv lv slab_thickness vacuum_thickness [asis] : step ga aru tame no model wo tsukuru (no edge)

--species --names A B C ... : Genso no namae wo sitei (kaku gensi, 9+ gyou me)
--species6 --names6 A B C ... : Genso no namae wo sitei (6 gyo me)
--noelement : Genso no namae wo 9 gyou me ikou hyouji sinai
--addnumber : Gensi tousi bangou wo 9 gyou me ikou hyouji suru

--incell : Subete no zahyou wo 0 to 1 no aida ni suru
--frac -f -d --fract --fractional --direct : Direct kara Cartesian ni suru
--cart -c : Cartesian kara Direct ni suru
--shift A B C : site wo zenbu (A,B,C) dake zurasu
--random [R] : gensi zahyouwo zenbu R ika ugokasu
--swap_atom A B : A-banme gensi to B-banme gensi no ichi dake wo koukansuru (genso ha koukan sinai)
--swap_species A B : A-banme genso to B-banme genso wo koukan suru
--similar POSCAR1 : POSCAR1 ni chikai zahyo no POSCAR wo kaesu
--nebimages POSCAR1 N : NEB you ni 00 kara N+1 made POSCAR image wo kaesu

--niggli:  POSCAR no niggli reduced cell wo tsukuru (POSCAR ha primitive to site atukau = primitive wo sagasanai!)
--reciniggli : POSCAR no gyaku kousi ga niggli reduced cell ni naru youni suru
--maxortho : POSCAR no c jiku wo a, b jiku ni taisi maximally orthogonal ni suru
--change_slab_vacuum_thickness thickness : slab model no sinkuu sou no atsusa wo kaeru (legacy option)
--change_slab_vacuum_thickness_ortho thickness : slab model no sinkuu sou no atsusa wo kaeru

--unique [tolerance] : POSCAR no tyouhuku suru saito wo hitotu ni suru
--makeslab minz maxz shiftz tolerance : z zahyou ga minz - maxz no sla wo tsukuru,	hituyou nara shiftz dake zurasu
--rm_atom A B C..D .. : gensi A B C..D .. wo kesu, gensi ha ikutu sitei sitemo ii, gensi ha 1 kara kazoeru
--rm_species A B C..D .. : gensi A B C..D .. wo kesu, gensi ha ikutu sitei sitemo ii, genso ha 1 kara kazoeru

--supercell A B C : POSCAR no supercell wo tukuru
--supercell A1 A2 A3 B1 B2 B3 C1 C2 C3 : POSCAR no supercell wo tukuru
--supercell [file] : POSCAR no supercell wo tukuru (matrix wo file ni ireru)
--add_atom x y z (namae) : POSCAR no saigo ni ichi (x,y,z) no gensi wo tasu, genso no namae ha (namae), namae ha shouryaku kanou	
--add_atom_exist x y z namae N : N ban me no genso no usiro ni (x,y,z) no gensi wo tasu,  namae ha shouryaku dekinai	
--recover_slab [isometry] : isometry de nai genshi wo oginau	
--uspex_slab h k l thickness : uspex hyoumen keisan you no model wo tsukuru
--addimages minx maxx miny maxy minz maxz : xyz muke ni gensi wo tasu

--cif : POSCAR wo P1 space group CIF ni suru, gensomei wa 6 gyome de uwagaki sinai!
--cif6 : POSCAR wo P1 space group CIF ni suru, gensomei wa 6 gyome de uwagaki!
--xyz : POSCAR wo xyz ni suru, gensomei wa 6 gyome de uwagaki sinai!
--xyz6 : POSCAR wo xyz ni suru, gensomei wa 6 gyome de uwagaki!

--cc --crystallographic_conventional --conventional --bposcar : Kessyougaku  *conventional* wo tsukuru
--cp --crystallographic_primtive : Kessyougaku  *primitive* wo tsukuru
--hinuma_reduced : Hinuma reduced triclinic cell wo tsukuru
--sc --standard_conventional : Standard *conventional* wo tsukuru
--sp --standard_primtive : Standard *primitive* wo tsukuru
--bzpoint [filename] : Brillouin zone no ten (kessyougaku ver) wo KPOINTS teki ni hyouji
--bzpoint_notrs [filename] : Brillouin zone no ten (no time reversal symmetry) wo KPOINTS teki ni hyouji
--scbzpoint [filename] : Brillouin zone no ten (SC ver) wo KPOINTS teki ni hyouji
--scpointkpf [filename] : Brillouin zone no ten (kessyougaku ver) wo kpf teki ni hyouji
--bzpointkpf [filename] : Brillouin zone no ten wo kpf teki ni hyouji
--bzpointkpf_notrs [filename] : Brillouin zone no ten (no time reversal symmetry) wo kpf teki ni hyouji
--scbzpointkpf [filename] : Brillouin zone no ten (SC ver) wo kpf teki ni hyouji
--kpath [interval] [filename] : suishyou kpath (kessyougaku ver) wo KPOINTS teki ni hyouji
--kpath_notrs [interval] [filename] : suishyou kpath (no time reversal symmetry) wo KPOINTS huu ni hyouji
--sckpath [interval] [filename] : suishyou kpath (SC ver) wo KPOINTS teki ni hyouji
--kpathkpf [filename] : suishyou kpath (kessyougaku ver) wo kpf teki ni hyouji
--kpathkpf_notrs [filename] : suishyou kpath (no time reversal symmetry) wo kpf teki ni hyouji
--sckpathkpf [filename] : suishyou kpath (SC ver) wo kpf teki ni hyouji
--kpathphonopy [filename] : suishyou kpath (kessyougaku ver) wo phonopy teki ni hyouji
--sckpathphonopy [filename] : suishyou kpath (SC ver) wo phonopy teki ni hyouji
--kpath_band [filename] : interval 0.025 kpath (kessyougaku ver) wo hyouji
--kpath_band_notrs [filename] : interval 0.025 kpath (no time reversal symmetry) wo hyouji
--sckpath_band [filename] : interval 0.025 kpath (SC ver) wo hyouji
--kpath_mass [filename] : yuukou situryou you 0.002 kpath (kessyougaku ver) wo hyouji
--kpath_mass_notrs [filename] : yuukou situryou you 0.002 kpath (no time reversal symmetry) wo hyouji
--sckpath_mass [filename] : yuukou situryou you 0.002 kpath (SC ver) wo hyouji
EOM
	exit;
}


sub version{
	print ("$version \n");
	exit;
}



#read/write 

sub readPOSCAR{
# $_[0], default STDIN no POSCAR wo yomikomu
# $_[1]="vasp4" nara 6 gyome wo yomanai
# Data wo hensuu to site hozon
# $topline: 1 gyou me
# $scale: 2 gyou me
# @latvec: 3-5 gyou me
# @namespecies: 6 gyou me
# @numspecies: 7 gyou me
# $dc: 8 gyou me (D or C)
# @x: x zahyou
# @y: y zahyou
# @z: z zahyou
# @w: saito no namae
# Tsuika no data
# $num_atoms: gensi no kazu
	$topline=&rmvcrlf("STDIN");
	my @a=&splitall("STDIN");
	$scale=$a[0];
	@latvec=&split3("STDIN");
	@latvec=(@latvec,&split3("STDIN"));
	@latvec=(@latvec,&split3("STDIN"));
	@namespecies=&splitall("STDIN");
	@numspecies=&splitall("STDIN");
	@a=&splitall("STDIN");
	if ($a[0]=~/^[Ss]/){
#selective dynamics: mou ikkai
		@a=&splitall("STDIN");
	}
	if ($a[0]=~/^[Dd]/){
		$dc="D";
	} elsif ($a[0]=~/^[CcKk]/){
		$dc="C";
	} else {
		die ("POSCAR no 8 gyoume ga DdCcKkSs no dore demo nai! \n");
	}
	$num_atoms=&numatoms;
	for (my $i=0; $i<$num_atoms; $i++){
		my @b=&split4("STDIN");
		$x[$i]=$b[0];
		$y[$i]=$b[1];
		$z[$i]=$b[2];
		$w[$i]=$b[3];
		if ($dc eq "C"){
			$x[$i]*=$scale;
			$y[$i]*=$scale;
			$z[$i]*=$scale;
		}
	}
	@latvec=&product_vs(@latvec,$scale);
	$scale=1;
}


sub clat{
#hikisuu no a b c alpha beta gamma kara lattice vector wo syuturyoku
#input: a b c alpha beta gamma (optional s)
#saigo ni "s" ga aru baai, hyouji sinai de atai wo kaesu (suppress)
	die ("--clat de hikisuu ga okasii desu: @_ \n") if ($_[5] eq ""); 
	my ($a, $b, $c, $al, $be, $ga)=@_[0..5];
	my @a=(&cosd($al),&cosd($be),&cosd($ga));
	my $vol2=1-&cosd($al)*&cosd($al)-&cosd($be)*&cosd($be)-&cosd($ga)*&cosd($ga);
	$vol2+=2*&cosd($al)*&cosd($be)*&cosd($ga);
	my $vol=sqrt($vol2);
	if ($_[6] ne "s"){
		&printx($a);
		&printx(0);
		&printx(0);
		print ("\n");
		&printx($b*&cosd($ga));
		&printx($b*&sind($ga));
		&printx(0);
		print ("\n");
		&printx($c*&cosd($be));
		&printx($c*(&cosd($al)-&cosd($be)*&cosd($ga))/&sind($ga));
		&printx($c*$vol/&sind($ga));
		print ("\n");
		exit;
	} else {
		return ($a,$b*&cosd($ga),$b*&sind($ga),$c*&cosd($be), $c*(&cosd($al)-&cosd($be)*&cosd($ga))/&sind($ga), $c*$vol/&sind($ga));
	}
}

sub writeCIF{
	&c2d if ($dc eq "C");
	print ("data_${topline}\n");
	print ("_pd_phase_name ${topline}\n");
	my @abcreal=&get_abc(@latvec);
	printf ("_cell_length_a  %2.10f\n",$abcreal[0]);
	printf ("_cell_length_b  %2.10f\n",$abcreal[1]);
	printf ("_cell_length_c  %2.10f\n",$abcreal[2]);
	printf ("_cell_angle_alpha  %2.10f\n",$abcreal[3]);
	printf ("_cell_angle_beta  %2.10f\n",$abcreal[4]);
	printf ("_cell_angle_gamma  %2.10f\n",$abcreal[5]);
	print ("_symmetry_space_group_name_H-M  'P1'\n");
	print ("_symmetry_Int_Tables_Number  1\n");
	print ("loop_\n");
	print ("_symmetry_equiv_pos_site_id\n");
	print ("_symmetry_equiv_pos_as_xyz_\n");
	print ("  1  x,y,z\n");
	print ("loop_\n");
	print (" _atom_site_label\n");
	print (" _atom_site_occupancy\n");
	print (" _atom_site_fract_x\n");
	print (" _atom_site_fract_y\n");
	print (" _atom_site_fract_z\n");
	print (" _atom_site_thermal_displace_type\n");
	print (" _atom_site_B_iso_or_equiv\n");
	print (" _atom_site_type_symbol\n");
	my $count=1;
	for (my $i=0; $i<=$#numspecies; $i++){
		my $name6=$namespecies[$i];
		$name6 = substr($name6, 0, 2) if (length($name6) > 2);
		$name6 =~ s/_//g;
		for (my $j=0; $j<$numspecies[$i]; $j++){
			my $namew=$w[$count-1];
			$namew = substr($namew, 0, 2) if (length($namew) > 2);
			if ($ARGV[0] eq "--cif6"){
				print("${name6}${count} 1.0 ");
			} else {
				print("${namew}${count} 1.0 ");
			}
			printf ("%2.10f ",$x[$count-1]);
			printf ("%2.10f ",$y[$count-1]);
			printf ("%2.10f ",$z[$count-1]);
			print(" Biso 1.0 ");
			if ($ARGV[0] eq "--cif6"){
				print("${name6}\n");
			} else {
				print("${namew}\n");
			}
			$count++;
		}
	}
	exit;
}

sub writexyz{
	&d2c if ($dc eq "D");
	print ("${num_atoms}\n");
	print ("data_${topline}\n");
	my $count=1;
	for (my $i=0; $i<=$#numspecies; $i++){
		my $name6=$namespecies[$i];
		$name6 = substr($name6, 0, 2) if (length($name6) > 2);
		for (my $j=0; $j<$numspecies[$i]; $j++){
			my $namew=$w[$count-1];
			$namew = substr($namew, 0, 2) if (length($namew) > 2);
			if ($ARGV[0] eq "--xyz6"){
				print("${name6} ");
			} else {
				print("${namew} ");
			}
			printf ("%2.10f ",$x[$count-1]);
			printf ("%2.10f ",$y[$count-1]);
			printf ("%2.10f \n",$z[$count-1]);
			$count++;
		}
	}
	exit;
}


sub writePOSCAR{
# POSCAR wo hyouji
# Ika no hensuu wo siyou
# $topline: 1 gyou me
# $scale: 2 gyou me
# @latvec: 3-5 gyou me
# @namespecies: 6 gyou me
# @numspecies: 7 gyou me
# $dc: 8 gyou me (D or C)
# @x: x zahyou
# @y: y zahyou
# @z: z zahyou
# @w: saito no namae
# $_: syuturyoku file: aru baai exit shinai
	my $handle="STDOUT";
	if ($_[0] ne ""){
		open OUT, "> $_[0]";
		$handle="OUT";
	}
	print $handle ("$topline\n");
	print $handle ("1.0\n");
	&printx($latvec[0],$handle);
	&printx($latvec[1],$handle);
	&printx($latvec[2],$handle);
	print $handle ("\n");
	&printx($latvec[3],$handle);
	&printx($latvec[4],$handle);
	&printx($latvec[5],$handle);
	print $handle ("\n");
	&printx($latvec[6],$handle);
	&printx($latvec[7],$handle);
	&printx($latvec[8],$handle);
	print $handle ("\n");
	print $handle ("@namespecies \n");
	print $handle ("@numspecies \n");
	if ($dc eq "D"){
		print $handle ("Direct(");
	} else {
		print $handle ("Cartesian(");
	}
	print $handle ("${num_atoms})\n");
	for (my $i=0; $i<$num_atoms; $i++){
		&printx($x[$i],$handle);
		&printx($y[$i],$handle);
		&printx($z[$i],$handle);
		if ($w[$i] eq ""){
			my @a=&nameatom($i+1);
			$w[$i]=$a[0];
		}
		print $handle ("  $w[$i] \n");
	}
	close OUT if ($_[0] ne "");
	exit if ($_[0] eq "");
}



sub writePOSCAR_selective{
# POSCAR wo hyouji
# Selective
# Ika no hensuu wo siyou
# $topline: 1 gyou me
# $scale: 2 gyou me
# @latvec: 3-5 gyou me
# @namespecies: 6 gyou me
# @numspecies: 7 gyou me
# $dc: 8 gyou me (D or C)
# @x: x zahyou
# @y: y zahyou
# @z: z zahyou
# @w: saito no namae
# $_: syuturyoku file: aru baai exit shinai
	print STDOUT ("$topline\n");
	print STDOUT ("1.0\n");
	&printx($latvec[0]);
	&printx($latvec[1]);
	&printx($latvec[2]);
	print STDOUT ("\n");
	&printx($latvec[3]);
	&printx($latvec[4]);
	&printx($latvec[5]);
	print STDOUT ("\n");
	&printx($latvec[6]);
	&printx($latvec[7]);
	&printx($latvec[8]);
	print STDOUT ("\n");
	print STDOUT ("@namespecies \n@numspecies \nSelective dynamics\n");
	if ($dc eq "D"){
		print STDOUT ("Direct(");
	} else {
		print STDOUT  ("Cartesian(");
	}
	print STDOUT ("${num_atoms})\n");
	for (my $i=0; $i<$num_atoms; $i++){
		&printx($x[$i]);
		&printx($y[$i]);
		&printx($z[$i]);
		if ($w[$i] eq ""){
			my @a=&nameatom($i+1);
			$w[$i]=$a[0];
		}
		print STDOUT (" $_[0] $_[0] $_[0] $w[$i] \n");
	}
	exit;
}



#4.14 keisiki de print
#hikisuu=(output, optional file handle)
sub printx{
	my $handle="STDOUT";
	$handle=$_[1] if ($_[1] ne "");
	my $a=$_[0];
	if ((int($a) eq 0) && ($a < 0)){
		printf $handle ("  %0.14f",$a);
	} else {
		printf $handle (" % 3d",$a);
		my $b=$a-int($a);
		$b*=-1 if ($b < 0);
		$b*=1E14;
		printf $handle (".%014d",$b);
	}
}

#file kara hikisuu kist wo nyuusyu
sub args_from_file{
	die ("file $_[0] ga arimasen") if (! -e $_[0]);
	open FILE, $_[0];
	my @out;
	while (my $a=<FILE>){
		@out=(@out, &splitall($a));
	}
	close FILE;
	return (@out);
}


#Individual subroutines
#Basic data


#Cartesian to Direct
sub c2d{
	for (my $i=0; $i<$num_atoms; $i++){
		my @invlatvec=&invmatrix3(@latvec);
		($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i],@invlatvec);
	}
}

#Direct to Cartesian 
sub d2c{
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i],@latvec);
	}
}

#Gensi suu
sub numatoms{
	my $a;
	for (my $i=0;$i<=$#numspecies;$i++){
		$a+=$numspecies[$i];
	}
	return ($a);
}

#ataerareta gensi no genso mei to bango wo kaesu
#saisyo no gensi ha 0
#saisyo no genso ha 1 (not 0)

sub nameatom{
#gensi no bangou
	my $num=$_[0];
	if ($num <= 0){
		die ("nameatom de gensi no bango $_[0] ga okasii ( <= 0 ) \n"); 
	}
	for (my $i=0; $i<=$#numspecies;$i++){
		$num-=$numspecies[$i];
		if ($num <= 0){
			return ($namespecies[$i],$i+1);
		}
	}
	die ("nameatom de gensi no bango $_[0] ga okasii (mosikasite jissai no gensi suu $num_atoms yori ookii? \n");
}

#gyakukousi vector
#hikisuu wo kousi teisuu to suru
#default: @latvec
sub reci_latvec{
	my @a=@_;
	@a=@latvec if ($_[0] eq "");
	my @reci1=(&det2(@a,4,5,7,8),&det2(@a,5,3,8,6),&det2(@a,3,4,6,7));
	my @reci2=(&det2(@a,7,8,1,2),&det2(@a,8,6,2,0),&det2(@a,6,7,0,1));
	my @reci3=(&det2(@a,1,2,4,5),&det2(@a,2,0,5,3),&det2(@a,0,1,3,4));
	my $det=&det3(@a);
	@reci1=&product_vs(@reci1,6.283185307179586/$det);
	@reci2=&product_vs(@reci2,6.283185307179586/$det);
	@reci3=&product_vs(@reci3,6.283185307179586/$det);
	for (my $i=0; $i<=2; $i++){
		$reci1[$i]=&precise($reci1[$i],0);
		$reci2[$i]=&precise($reci2[$i],0);
		$reci3[$i]=&precise($reci3[$i],0);
	}
	my @reci=(@reci1, @reci2, @reci3);
	return (@reci);
}

#kousi teisuu wo hikisuu kara nyuusyu
#a b c alpha beta gamma wo kaesu (degrees)
#scale wa kangaenai
sub get_abc{
	my $a=&norm(@_[0..2]);
	my $b=&norm(@_[3..5]);
	my $c=&norm(@_[6..8]);
	my $al=&acos(&product_vv(@_[3..5],@_[6..8])/$b/$c);
	my $be=&acos(&product_vv(@_[6..8],@_[0..2])/$c/$a);
	my $ga=&acos(&product_vv(@_[0..2],@_[3..5])/$a/$b);
	$b=&precise($b,$a);
	$c=&precise($c,$a);
	$c=&precise($c,$b);
	$al=&precise($al,60);
	$be=&precise($be,60);
	$ga=&precise($ga,60);
	$al=&precise($al,90);
	$be=&precise($be,90);
	$ga=&precise($ga,90);
	$al=&precise($al,120);
	$be=&precise($be,120);
	$ga=&precise($ga,120);
	$be=&precise($be,$al);
	$ga=&precise($ga,$al);
	$ga=&precise($ga,$be);
	return ($a, $b, $c, $al, $be, $ga);
}


#tyouhuku site wo matomeru
#Direct ni suru
sub unique{
	my $tolerance=1E-8;
	$tolerance=$_[0] if ($_[0] ne "");
	my (@x_new, @y_new, @z_new,@w_new);
#genso goto ni shori
	for (my $i=0;$i<=$#numspecies;$i++){
#genso goto ni @x_work etc wo tukuru 
		my @x_work=splice(@x,0,$numspecies[$i]);
		my @y_work=splice(@y,0,$numspecies[$i]);
		my @z_work=splice(@z,0,$numspecies[$i]);
		my @w_work=splice(@w,0,$numspecies[$i]);
#tyouhuku wo korosu
		my @work_new=&unique_species($tolerance,$numspecies[$i],@x_work,@y_work,@z_work,@w_work);
		$numspecies[$i]=shift @work_new;
		my @x_work_new=splice(@work_new,0,$numspecies[$i]);
		my @y_work_new=splice(@work_new,0,$numspecies[$i]);
		my @z_work_new=splice(@work_new,0,$numspecies[$i]);
		my @w_work_new=splice(@work_new,0,$numspecies[$i]);
#atarasii hairetu
		@x_new=(@x_new,@x_work_new);
		@y_new=(@y_new,@y_work_new);
		@z_new=(@z_new,@z_work_new);
		@w_new=(@w_new,@w_work_new);
	}
#gensisuu wo update
	$num_atoms=&numatoms;
#moto no hairetu wo update
	@x=@x_new;
	@y=@y_new;
	@z=@z_new;
	@w=@w_new;
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
}

#tyouhuku site wo matomeru
#Nyuuryoku ha hairetu (tolerance, gensisuu, @x, @y, @z, @w) (not unique)
#Hairetu (gensisuu, @x, @y, @z, @w) (unique) wo kaesu
#Note: hitotu no genso dake wo kangaeru!
sub unique_species{
#hensuu
	my @in=@_;
#moto no gensi
	my $tolerance=shift @in;
	my $natoms_in=shift @in;
	my @x_in=splice(@in,0,$natoms_in);
	my @y_in=splice(@in,0,$natoms_in);
	my @z_in=splice(@in,0,$natoms_in);
	my @w_in=splice(@in,0,$natoms_in);
#atarasii no gensi
	my $natoms_out=0;
	my (@x_out, @y_out, @z_out,@w_out);
#tyouhuku kakunin
	for (my $i=0;$i<$natoms_in;$i++){
		my $duplicate_flag=0;
		for (my $j=0; $j<$natoms_out; $j++){
#tyouhuku ga aruto $duplicate_flag ga 0 de nakunaru
			$duplicate_flag+=&site_is_same($tolerance,$x_in[$i],$y_in[$i],$z_in[$i],$x_out[$j],$y_out[$j],$z_out[$j]);
		}
		if ($duplicate_flag == 0){
#tyouhuku nasi
			$x_out[$natoms_out]=$x_in[$i];
			$y_out[$natoms_out]=$y_in[$i];
			$z_out[$natoms_out]=$z_in[$i];
			$w_out[$natoms_out]=$w_in[$i];
			$natoms_out++;
		}
	}
	return ($natoms_out,@x_out,@y_out,@z_out,@w_out);
}


#tyokuhoutai supercell wo tukuru
#kyoukai wo (xmin,xmax,ymin,ymax,zmin,zmax) no you ni ireru
#Direct ni suru koto!
#gensi ha kezura nai
sub populate_block{
	my @a=@_;
#keisuu
	my $minx=&floor(shift @a);
	my $maxx=&ceiling(shift @a);
	my $miny=&floor(shift @a);
	my $maxy=&ceiling(shift @a);
	my $minz=&floor(shift @a);
	my $maxz=&ceiling(shift @a);
	my $delx=$maxx-$minx;
	my $dely=$maxy-$miny;
	my $delz=$maxz-$minz;
	my $scale=$delx*$dely*$delz;
	my (@xnew, @ynew, @znew,@wnew);
#kousi teisuu
	@latvec[0..2]=&product_vs(@latvec[0..2],$delx);
	@latvec[3..5]=&product_vs(@latvec[3..5],$dely);
	@latvec[6..8]=&product_vs(@latvec[6..8],$delz);
	#gensi wo update
	for (my $i=0;$i<$num_atoms;$i++){
		for (my $j=$minx;$j<$maxx;$j++){
			for (my $k=$miny;$k<$maxy;$k++){
				for (my $l=$minz;$l<$maxz;$l++){
					push (@xnew,($x[$i]+$j)/$delx);
					push (@ynew,($y[$i]+$k)/$dely);
					push (@znew,($z[$i]+$l)/$delz);
					push (@wnew,$w[$i]);
				}
			}
		}
	}
	@x=@xnew;
	@y=@ynew;
	@z=@znew;
	@w=@wnew;
#gensisuu wo update
	$num_atoms=0;
	for (my $i=0;$i<=$#numspecies;$i++){
		$numspecies[$i]*=$scale;
		$num_atoms+=$numspecies[$i];
	}
}

# --data
# --dataniggli demo tsukau
sub data{
	print ("$version");
	print ("REAL LATTICE\n");
	print (" Real space a b c alpha beta gamma: ");
	my @abcreal=&get_abc(@latvec);
	if ($_[0] eq "full"){
		printf ("%2.10f %2.10f %2.10f %3.10f %3.10f %3.10f\n",@abcreal[0..5]); 
	} else {
		printf ("%2.10f %2.10f %2.10f %3.3f %3.3f %3.3f\n",@abcreal[0..5]); 
	}
	print (" Real space Volume: ");
	my @ab=&cross_vv(@latvec[0..5]);
	my $vol=&product_vv(@ab,@latvec[6..8]);
	if ($_[0] eq "full"){
		printf ("%4.10f\n",$vol); 
	} else {
		printf ("%4.4f\n",$vol); 
	}
	print (" Real space c/a = ");
	if ($_[0] eq "full"){
		printf ("%1.10f ",$abcreal[2]/$abcreal[0]); 
	} else {
		printf ("%1.4f ",$abcreal[2]/$abcreal[0]); 
	}
	print (" |axb| = ");
	my $abnorm=&norm(@ab);
	if ($_[0] eq "full"){
		printf ("%2.10f ",$abnorm); 
	} else {
		printf ("%2.4f ",$abnorm); 
	}
	print (" vol/|axb| = ");
	if ($_[0] eq "full"){
		printf ("%2.10f\n",$vol/$abnorm); 
	} else {
		printf ("%2.4f\n",$vol/$abnorm); 
	}
	print ("RECIPROCAL LATTICE\n");
	print (" Reciprocal space lattice:\n");
	my @recilatvec=&reci_latvec(@latvec);
	printf ("   %2.10f  %2.10f  %2.10f\n",@recilatvec[0..2]); 
	printf ("   %2.10f  %2.10f  %2.10f\n",@recilatvec[3..5]); 
	printf ("   %2.10f  %2.10f  %2.10f\n",@recilatvec[6..8]); 
	print (" Reciprocal space a b c alpha beta gamma: ");
	my @abcreci=&get_abc(@recilatvec);
	if ($_[0] eq "full"){
		printf ("%2.10f %2.10f %2.10f %3.10f %3.10f %3.10f\n",@abcreci); 	} else {
		printf ("%2.10f %2.10f %2.10f %3.3f %3.3f %3.3f\n",@abcreci); 
	}
	print (" Reciprocal space Volume: ");
	if ($_[0] eq "full"){
		printf ("%4.10f\n",&det3(@recilatvec)); 
 	} else {
		printf ("%4.4f\n",&det3(@recilatvec)); 
	}
	exit;
}

#gensi ga nai genso wo kesu
sub remove_zero_atom_species{
	for (my $i=$#numspecies; $i>=0;$i--){
		if ($numspecies[$i] == 0){
			splice (@namespecies, $i, 1);
			splice (@numspecies, $i, 1);
		}
	}
}

sub distance{
#kyori wo hyouji
#+/- $_[0] no supercell wo kentou
#$_[1] ga "default": $_[2] ika no subete no kyori wo hyouji
#$_[1] ga "min" : saisyou kyori wo hyoji
#$_[1] ga "coordination" : $_[2], $_[3] wo moto ni coordination wo hyouji

#kyori wo siraberu
#hanni
	my $max_dist;
	my $cutoff1=$_[2];
	my $cutoff2=$_[3];
	if ($_[1] eq "default"){
		$max_dist=$_[2];
		$max_dist=4 if ($max_dist eq "");
	} elsif ($_[1] eq "min"){
		$max_dist=9999;
	} elsif ($_[1] eq "coordination"){
		$max_dist=4;
		$cutoff1=1.1 if ($cutoff1 eq "");
		$cutoff2=1.18 if ($cutoff2 eq "");
	} elsif ($_[1] eq "list") {
		$max_dist=$_[2];
	}
	my @latvec_bak=@latvec;
#list ni tsukau
	my @argument=@_;
	shift @argument;
	shift @argument;
	shift @argument;
#header
	if ($_[1] eq "default"){
		print ("$version");
		print ("Gensi no sitei ha\n");
		print ("[Gensi no bangou] [namae] [n1 n2 n3 (unit cell)]\n");
		print ("kyori ha $arg2[0] no hanni de sagasu\n\n");
		print ("******************** KYORI ********************\n");
	}
#hikae
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @w_orig=@w;
	my @latvec_orig=@latvec;
#dokomade siraberuka
	my $scalemin=-1*($_[0]);
	my $scalemax=($_[0])+1;
#ikkai D ni suru
	&c2d if ($dc eq "C");
#hikae wo toru
	my @x_orig_d=@x;
	my @y_orig_d=@y;
	my @z_orig_d=@z;
#scale
	my (@x_small, @y_small, @z_small);
	for (my $i=0; $i<$num_atoms; $i++){
		$x_small[$i]=$x[$i]/($scalemax-$scalemin);
		$y_small[$i]=$y[$i]/($scalemax-$scalemin);
		$z_small[$i]=$z[$i]/($scalemax-$scalemin);
	}
	&populate_block($scalemin,$scalemax,$scalemin,$scalemax,$scalemin,$scalemax);
#kono dankai de @x_small toka ha D de moto no gensi saito ga umatteru
#@x toka ha D de image saito mo umatteru
#@x toka wo C ni suru
	&d2c;
	my @x_large_c=@x;
	my @y_large_c=@y;
	my @z_large_c=@z;
	my @w_large=@w;
	my @w_near; 
#@x_small toka wo @x toka ni site C ni suru
	@x=@x_small;
	@y=@y_small;
	@z=@z_small;
	@w=@w_orig;
	&d2c;
#kokode @x toka ha moto no gensi saito dake umatteru
#@x_large toka ha image mo umatteru
#ryouhou C
#saito goto ni kyori wo dasite sort
	for (my $i=0; $i<=$#w_orig; $i++){
		if ($_[1] eq "default"){
			printf ("%-5s%-7s%1.4f  %1.4f  %1.4f\n", $i+1,$w_orig[$i],$x_orig[$i],$y_orig[$i],$z_orig[$i]) ;
		} elsif (($_[1] eq "list") && (grep { $_ == ($i+1) } @argument )){
			printf ("%-5s%-7s%1.4f  %1.4f  %1.4f\n", $i+1,$w_orig[$i],$x_orig[$i],$y_orig[$i],$z_orig[$i]) ;
		}
		my @x_large_copy=@x_large_c;
		my @y_large_copy=@y_large_c;
		my @z_large_copy=@z_large_c;
		my (@x_near, @y_near, @z_near, @w_near, @n_near, @dist_near);
		for (my $m=0;$m<=$#w_orig;$m++){
			for (my $j=$scalemin;$j<$scalemax;$j++){
				for (my $k=$scalemin;$k<$scalemax;$k++){
					for (my $l=$scalemin;$l<$scalemax;$l++){
						my @shifted=(shift(@x_large_copy),shift(@y_large_copy),shift(@z_large_copy));
						my @vector=&diff_vv(@shifted,$x[$i],$y[$i],$z[$i]);
						my $dist=&norm(@vector);
						$max_dist=$dist if (($dist<$max_dist) && ($dist > 0.00001) && ($_[1] eq "min"));
						if ($dist < $max_dist){
							if (($dist > 0.00001) || ($i != $m)){
								push (@x_near,$j);
								push (@y_near,$k);
								push (@z_near,$l);
								push (@n_near,$m);
								push (@w_near,$w_orig[$m]);
								push (@dist_near,$dist);

							}
						}
					}
				}
			}
		}
		my @order_key=&order_key_number(@dist_near);
		@x_near=@x_near[@order_key];
 		@y_near=@y_near[@order_key];
		@z_near=@z_near[@order_key];
		@n_near=@n_near[@order_key];
		@w_near=@w_near[@order_key];
		@dist_near=@dist_near[@order_key];
#print normal
		if ($_[1] eq "default"){
			for (my $n=0; $n <= $#dist_near; $n++){
				printf ("      %-5s%-6s%2s%4s%4s    %2.4f\n",$n_near[$n]+1,$w_near[$n],$x_near[$n],$y_near[$n],$z_near[$n],$dist_near[$n]);
			}
		} elsif (($_[1] eq "list") && (grep { $_ == ($i+1) } @argument )){
			for (my $n=0; $n <= $#dist_near; $n++){
				printf ("      %-5s%-6s%2s%4s%4s    %2.4f\n",$n_near[$n]+1,$w_near[$n],$x_near[$n],$y_near[$n],$z_near[$n],$dist_near[$n]);
			}
		} elsif ($_[1] eq "coordination"){
			printf ("%-5s%-7s ", $i+1,$w_orig[$i]);
			my $coordination=0;
			for (my $n=0; $n <= $#dist_near; $n++){
				$coordination+=0.3 if ($dist_near[$n] < ($dist_near[0]*$cutoff1));
				$coordination+=0.7 if ($dist_near[$n] < ($dist_near[0]*$cutoff2));
			}
			print ("$coordination ");
			if (($coordination == 2) || ($coordination == 3)){
				my @neighborlist;
				for (my $n=0; $n <= $#dist_near; $n++){
					if ($dist_near[$n] < ($dist_near[0]*$cutoff1)){
						my @neighbor=&product_vm($x_near[$n],$y_near[$n],$z_near[$n],@latvec_bak);
						push (@neighborlist,&sum_vv(@neighbor,$x[$n_near[$n]],$y[$n_near[$n]],$z[$n_near[$n]]));
					} 
				}
				printf ("%3.2f ",&angle_vv($x[$i],$y[$i],$z[$i],@neighborlist[0..5]));
				if ($coordination == 3){
					printf ("%3.2f ",&angle_vv($x[$i],$y[$i],$z[$i],@neighborlist[0..2],@neighborlist[6..8]));
					printf ("%3.2f ",&angle_vv($x[$i],$y[$i],$z[$i],@neighborlist[3..8]));
				}
			}
			print ("\n");
		}
	}
	print ("$max_dist \n") if ($_[1] eq "min");
	exit;
}


sub distance_sphere{
#kyori $_[0] nai no gensi ichi wo unit  ni shaei
	my $max_dist=$_[0];
#hikae
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @latvec_orig=@latvec;

#scale
	my $abnorm=&norm(&cross_vv(@latvec[0..5]));
	my $acnorm=&norm(&cross_vv(@latvec[0..2],@latvec[6..8]));
	my $bcnorm=&norm(&cross_vv(@latvec[3..8]));
	my $vol=&det3(@latvec[0..8]);
	my $amin=-&ceiling($max_dist/$vol*$bcnorm+1);
	my $bmin=-&ceiling($max_dist/$vol*$acnorm+1);
	my $cmin=-&ceiling($max_dist/$vol*$abnorm+1);
	my $amax=1-$amin;
	my $bmax=1-$bmin;
	my $cmax=1-$cmin;
#ikkai D ni suru
	&c2d if ($dc eq "C");
#hikae wo toru
	my @x_orig_d=@x;
	my @y_orig_d=@y;
	my @z_orig_d=@z;
#small wo hozon (moto no cell wo supercell no zahyou)
	my (@x_small, @y_small, @z_small);
	for (my $i=0; $i<$num_atoms; $i++){
		$x_small[$i]=$x[$i]/($amax-$amin);
		$y_small[$i]=$y[$i]/($bmax-$bmin);
		$z_small[$i]=$z[$i]/($cmax-$cmin);
	}
	&populate_block($amin,$amax,$bmin,$bmax,$cmin,$cmax);

#kono dankai de @x_small toka ha D de moto no gensi saito ga umatteru
#@x toka ha D de image saito mo umatteru
#@x toka wo C ni suru
#large wo hozon (supercell wosupercell no zahyou)
	&d2c;
	my @x_large_c=@x;
	my @y_large_c=@y;
	my @z_large_c=@z;
#@x_small toka wo @x toka ni site C ni suru
	@x=@x_small;
	@y=@y_small;
	@z=@z_small;
	&d2c;
#kokode @x toka ha moto no gensi saito dake umatteru
#@x_large toka ha image mo umatteru
#ryouhou C
#saito goto ni kyori wo dasite sort
	for (my $i=0; $i<=$#x_orig;$i++){
		my $tt=$i+1;
		print ("atom $tt ");
		my @x_large_copy=@x_large_c;
		my @y_large_copy=@y_large_c;
		my @z_large_copy=@z_large_c;
		my (@xx,@yy,@zz);
		for (my $m=0;$m<=$#x_orig;$m++){
			for (my $j=$amin;$j<$amax;$j++){
				for (my $k=$bmin;$k<$bmax;$k++){
					for (my $l=$cmin;$l<$cmax;$l++){
						my @shifted=(shift(@x_large_copy),shift(@y_large_copy),shift(@z_large_copy));
						my @vector=&diff_vv(@shifted,$x[$i],$y[$i],$z[$i]);
#normalize
						my $dist=&norm(@vector);
						if (&within($dist,0.0001,$max_dist)){
							@vector_norm=&product_vs(@vector,1/$dist);

							push (@xx,$vector_norm[0]);
							push (@yy,$vector_norm[1]);
							push (@zz,$vector_norm[2]);
						}
					}
				}
			}
		}
#remove duplicate
		for (my $m=$#xx; $m>0;$m--){
			for (my $n=$m-1; $n>=0;$n--){
				my $flag=0;
				my $t=&norm(&diff_vv($xx[$m],$yy[$m],$zz[$m],$xx[$n],$yy[$n],$zz[$n]));
				if (($t < 1E-8)&&($flag==0)){
					$flag=1;
					$xxx=splice(@xx,$m,1);
					$yyy=splice(@yy,$m,1);
					$zzz=splice(@zz,$m,1);
				}
			}
		}
# ga deta node qhull no jyunbi
		my @r0=(&part_sum_array($#xx,@xx),&part_sum_array($#yy,@yy),&part_sum_array($#zz,@zz));
		@r0=&product_vs(@r0,-1/($#xx+1));
		my $r0norm=&norm(@r0);
		if ($r0norm < 0.01){
			print ("in_bulk\n");
		} else {
			@r0=&product_vs(@r0,1/$r0norm);
			unshift(@xx,$r0[0]);
			unshift(@yy,$r0[1]);
			unshift(@zz,$r0[2]);
			open QHULL_IN, ">tsubo_temp_qhull_in";
			my $t=$#xx+1;
			print QHULL_IN ("3\n$t\n");
			for (my $m=0; $m<=$#xx; $m++){
				print QHULL_IN ("$xx[$m] $yy[$m] $zz[$m]\n");
			}
			close QHULL_IN;
			`cat tsubo_temp_qhull_in | qconvex o >tsubo_temp_qhull_out`;
			open QHULL_OUT, "tsubo_temp_qhull_out";
			my $a=<QHULL_OUT>;
			$a=<QHULL_OUT>;
			my @b=&splitall($a);
			my $sum=0;
			for (my $m=1; $m<=$b[0]; $m++){
				$a=<QHULL_OUT>;
			}
			for (my $m=1; $m<=$b[1]; $m++){
				$a=<QHULL_OUT>;
				my @c=&splitall($a);
				if (($c[4] ne "") && ($c[1]*$c[2]*$c[3]*$c[4] == 0)){
					die ("qhull de spherical triange igai wo tasu koto ni naru: gensi no zahyou wo zurasu hituyou ari!? \n")
				}
				if (($c[1]*$c[2]*$c[3]) == 0){
					my @angle;
					$angle[0]=&product_vv($xx[$c[1]],$yy[$c[1]],$zz[$c[1]],$xx[$c[2]],$yy[$c[2]],$zz[$c[2]]);
					$angle[1]=&product_vv($xx[$c[1]],$yy[$c[1]],$zz[$c[1]],$xx[$c[3]],$yy[$c[3]],$zz[$c[3]]);
					$angle[2]=&product_vv($xx[$c[2]],$yy[$c[2]],$zz[$c[2]],$xx[$c[3]],$yy[$c[3]],$zz[$c[3]]);
					$sum+=&spherical_triangle_area(@angle);
				}
			}
			$sum/=(atan2(1,1)*4);
			print ("$sum \n");
		}
	}
	`rm tsubo_temp_qhull_in tsubo_temp_qhull_out`;
	exit;
}




# jiku wo kaiten
# yurusareru hikisuu ha
# abc bca cab -acb cb-a b-ac
# c-ba -bac ac-b ba-c a-cb -cba
# a-b-c -b-ca -ca-b -ab-c b-c-a -c-ab
# -a-bc -bc-a c-a-b -a-c-b -c-b-a -b-a-c
# mosikuha 1~24 

sub rotate_axis{
	my @a=@latvec[0..2];
	my @b=@latvec[3..5];
	my @c=@latvec[6..8];
	my @d=&product_vs(@latvec[0..2],-1);
	my @e=&product_vs(@latvec[3..5],-1);
	my @f=&product_vs(@latvec[6..8],-1);
	my @g=("hikisu_ga_chigau");
	@g=(@a,@b,@c) if (($_[0] eq "abc") || ($_[0]==1));
	@g=(@b,@c,@a) if (($_[0] eq "bca") || ($_[0]==2));
	@g=(@c,@a,@b) if (($_[0] eq "cab") || ($_[0]==3));
	@g=(@d,@c,@b) if (($_[0] eq "-acb") || ($_[0]==4));
	@g=(@c,@b,@d) if (($_[0] eq "cb-a") || ($_[0]==5));
	@g=(@b,@d,@c) if (($_[0] eq "b-ac") || ($_[0]==6));
	@g=(@c,@e,@a) if (($_[0] eq "c-ba") || ($_[0]==7));
	@g=(@e,@a,@c) if (($_[0] eq "-bac") || ($_[0]==8));
	@g=(@a,@c,@e) if (($_[0] eq "ac-b") || ($_[0]==9));
	@g=(@b,@a,@f) if (($_[0] eq "ba-c") || ($_[0]==10));
	@g=(@a,@f,@b) if (($_[0] eq "a-cb") || ($_[0]==11));
	@g=(@f,@b,@a) if (($_[0] eq "-cba") || ($_[0]==12));
	@g=(@a,@e,@f) if (($_[0] eq "a-b-c") || ($_[0]==13));
	@g=(@e,@f,@a) if (($_[0] eq "-b-ca") || ($_[0]==14));
	@g=(@f,@a,@e) if (($_[0] eq "-ca-b") || ($_[0]==15));
	@g=(@d,@b,@f) if (($_[0] eq "-ab-c") || ($_[0]==16));
	@g=(@b,@f,@d) if (($_[0] eq "b-c-a") || ($_[0]==17));
	@g=(@f,@d,@b) if (($_[0] eq "-c-ab") || ($_[0]==18));
	@g=(@d,@e,@c) if (($_[0] eq "-a-bc") || ($_[0]==19));
	@g=(@e,@c,@d) if (($_[0] eq "-bc-a") || ($_[0]==20));
	@g=(@c,@d,@e) if (($_[0] eq "c-a-b") || ($_[0]==21));
	@g=(@d,@f,@e) if (($_[0] eq "-a-c-b") || ($_[0]==22));
	@g=(@f,@e,@d) if (($_[0] eq "-c-b-a") || ($_[0]==23));
	@g=(@e,@d,@f) if (($_[0] eq "-b-a-c") || ($_[0]==24));
	die ("rotate_axis no hikisuu ga okasii \n") if ($g[0] eq "hikisu_ga_chigau");	@latvec=@g;
}


sub get_additional_coordinates{
#tuika POSCAR no zahyo wo nyushu
#hikisuu: POSCAR no ichi
#$_[1] eq "no_atom_check": gensi suu check wo hazusu

#(x1, y1, z1, x2, y2, z2...) wo hairetsu to shite kaesu
	die ("get_additional_coordinates: yomikomu file ga sonzai sinai \n") if (! -e $_[0]);
#kousi teisuu
	open POSCAR1,$_[0];
	my $a=<POSCAR1>;
	my $scale2=&rmvcrlf($a=<POSCAR1>);
	my @latvec2=&split3($a=<POSCAR1>);
	@latvec2=(@latvec2,&split3($a=<POSCAR1>));
	@latvec2=(@latvec2,&split3($a=<POSCAR1>));
	@latvec2=&product_vs(@latvec2, $scale2);
	$a=<POSCAR1>;
#gensisuu
	my @numspecies2=&splitall($a=<POSCAR1>);
	die ("species no kazu ga chigau: STDIN, $_[0] \n") if ($#numspecies ne $#numspecies2);
#debug 2680/4/12
	my @diffnumspecies=&diff_vv(@numspecies, @numspecies2);
	my $normdiffnumspecies=&norm(@diffnumspecies);
	die ("gensi suu ga chigau:STDIN, $_[0] \n") if (($normdiffnumspecies != 0) && ($_[1] ne "no_atom_check"));
	$a=<POSCAR1>;
	die ("$_[0] ga direct de nai") if ($a !~ /^[Dd]/);
#naibu zahyou
	my @list;
	for (my $i=0;$i<$num_atoms;$i++){
		$a=<POSCAR1>;
		my @b=&split3($a);
		push (@list, $b[0]);
		push (@list, $b[1]);
		push (@list, $b[2]);
	}
	close POSCAR1;
	return (@latvec2, @list);
}



#supercell wo tsukuru
sub get_supercell{
#matrix wo kakunin
	my @matrix=&get_transformation_matrix(@_);
#Determinant 0 wo hajiku
	my $det=abs(&det3(@matrix[0..8]));
	die ("@matrix no determinant ga 0 \n") if ($det == 0);

#Direct ni suru
	&c2d if ($dc eq "C");
#numatoms check
	my $num_atoms_orig=$num_atoms;
#tyuokuhoutai cell wo tukuru
	my @latvec_orig=@latvec;
	my $minx=&min(0,$matrix[0])+&min(0,$matrix[3])+&min(0,$matrix[6]);
	my $miny=&min(0,$matrix[1])+&min(0,$matrix[4])+&min(0,$matrix[7]);
	my $minz=&min(0,$matrix[2])+&min(0,$matrix[5])+&min(0,$matrix[8]);
	my $maxx=&max(0,$matrix[0])+&max(0,$matrix[3])+&max(0,$matrix[6]);
	my $maxy=&max(0,$matrix[1])+&max(0,$matrix[4])+&max(0,$matrix[7]);
	my $maxz=&max(0,$matrix[2])+&max(0,$matrix[5])+&max(0,$matrix[8]);
	&populate_block($minx,$maxx,$miny,$maxy,$minz,$maxz);
#Cart ni suru
	&d2c;
#Koushi teisuu wo tadasii mono ni kaeru
	@latvec=&product_mm(@matrix[0..8],@latvec_orig);
#Direct ni suru
	&c2d;
	&unique();
#numatoms check
	my $num_atoms_new=$num_atoms;
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i]+$matrix[9],$y[$i]+$matrix[10],$z[$i]+$matrix[11]);
	}
	die ("Supercell no gensisuu ga okasii \n") if (($num_atoms_new/$det) != $num_atoms_orig);
	&d2c if ($dc eq "C");
}



#niggli reduced cell no lattice vector wo motomeru
#algorithm in Acta Cryst (1976) A 32, 297 Krivy and Gruber
sub get_niggli_latvec{
	my $tolerance=1E-12;
	my $a=$_[0]*$_[0]+$_[1]*$_[1]+$_[2]*$_[2];
	my $b=$_[3]*$_[3]+$_[4]*$_[4]+$_[5]*$_[5];
	my $c=$_[6]*$_[6]+$_[7]*$_[7]+$_[8]*$_[8];
	my $p=2*($_[3]*$_[6]+$_[4]*$_[7]+$_[5]*$_[8]);
	my $q=2*($_[0]*$_[6]+$_[1]*$_[7]+$_[2]*$_[8]);
	my $r=2*($_[3]*$_[0]+$_[4]*$_[1]+$_[5]*$_[2]);
	my @vec1=@_[0..2];
	my @vec2=@_[3..5];
	my @vec3=@_[6..8];
	my $flag=1;
	my $back=0;
	while ($flag > 0){
		$a=0 if (abs($a) < $tolerance);
		$b=0 if (abs($b) < $tolerance);
		$c=0 if (abs($c) < $tolerance);
		$p=0 if (abs($p) < $tolerance);
		$q=0 if (abs($q) < $tolerance);
		$r=0 if (abs($r) < $tolerance);
		
		$flag=0;
#step 1
		$flag++ if (($a-$b) > $tolerance);
		$flag++ if ((abs($a-$b) < $tolerance) && ((abs($p)-abs($q)) > $tolerance));
		if ($flag > 0){
			($a,$b)=($b,$a);
			($p,$q)=($q,$p);
			my @w=@vec1;
			@vec1=@vec2;
			@vec2=@w;
		}
#step 2
		$flag=0;
		$flag++ if (($b-$c) > $tolerance);
		$flag++ if ((abs($b-$c) < $tolerance) && ((abs($q)-abs($r)) > $tolerance));
		if ($flag > 0){
			($b,$c)=($c,$b);
			($q,$r)=($r,$q);
			my @w=@vec2;
			@vec2=@vec3;
			@vec3=@w;
			next;
		}
#step 3,4
		if (($p*$q*$r) > $tolerance){
#step 3: flip a vector if necessary
			@vec1=&product_vs(@vec1,-1) if (($q<0) && ($r<0));
			@vec2=&product_vs(@vec2,-1) if (($r<0) && ($p<0));
			@vec3=&product_vs(@vec3,-1) if (($p<0) && ($q<0));
			($p,$q,$r)=(abs($p),abs($q), abs($r));
		} else {
#step 4: count number of negatives
			my $countn=0;
			$countn++ if ($p < -$tolerance);
			$countn++ if ($q < -$tolerance);
			$countn++ if ($r < -$tolerance);
			my $count0=0;
			$count0++ if ((abs($p)) < $tolerance);
			$count0++ if ((abs($q)) < $tolerance);
			$count0++ if ((abs($r)) < $tolerance);
			if ($countn==1){
#one negative: flip an axies
				if ($p < -$tolerance) {
					@vec1=&product_vs(@vec1,-1);
				} elsif ($q < -$tolerance) {
					@vec2=&product_vs(@vec2,-1);
				} else {
					@vec3=&product_vs(@vec3,-1);
				}
			} elsif ($countn==0){
#Zero negative: flip an axis if one or two positives
				if ($count0 == 2){
#one positive
					if ($p > $tolerance) {
						@vec2=&product_vs(@vec2,-1);
					} elsif ($q > $tolerance) {
						@vec3=&product_vs(@vec3,-1);
					} else {
						@vec1=&product_vs(@vec1,-1);
					}
#two positive
				 } elsif ($count0 == 1){
					if ($p < $tolerance) {
						@vec1=&product_vs(@vec1,-1);
					} elsif ($q < $tolerance) {
						@vec2=&product_vs(@vec2,-1);
					} else {
						@vec3=&product_vs(@vec3,-1);
					}
				}
			}
			($p,$q,$r)=(-1*abs($p),-1*abs($q), -1*abs($r));
		}
#step 5
		$flag=0;
		$flag++ if ((abs($p)-$b) > $tolerance);
		$flag++ if ((abs($p-$b) < $tolerance) && (($r-(2*$q)) > $tolerance));
		$flag++ if ((abs($p+$b) < $tolerance) && ($r<-$tolerance));
		if ($flag > 0){
			if ($p<0){
				@vec3=&sum_vv(@vec3,@vec2);
			} else {
				@vec3=&diff_vv(@vec3,@vec2);
			}
			my $sign=$p/abs($p);
			$c=$b+$c-$p*$sign;
			$q=$q-$r*$sign;
			$p=$p-2*$b*$sign;
			$back=5;
			next;
		}
#step 6
		$flag=0;
		$flag++ if ((abs($q)-$a) > $tolerance);
		$flag++ if ((abs($q-$a) < $tolerance) && (($r-(2*$p)) > $tolerance));
		$flag++ if ((abs($q+$a) < $tolerance) && ($r<-$tolerance));
		if ($flag > 0){
			if ($q<0){
				@vec3=&sum_vv(@vec3,@vec1);
			} else {
				@vec3=&diff_vv(@vec3,@vec1);
			}
			my $sign=$q/abs($q);
			$c=$a+$c-$q*$sign;
			$p=$p-$r*$sign;
			$q=$q-2*$a*$sign;
			$back=6;
			next;
		}
#step 7
		$flag=0;
		$flag++ if ((abs($r)-$a) > 0.0000091);
		$flag++ if ((abs($r-$a) < $tolerance) && (($q-(2*$r)) > $tolerance));
		$flag++ if ((abs($r+$a) < $tolerance) && ($q<-$tolerance));
		if ($flag > 0){
			if ($r<0){
				@vec2=&sum_vv(@vec2,@vec1);
			} else {
				@vec2=&diff_vv(@vec2,@vec1);
			}
			my $sign=$r/abs($r);
			$b=$a+$b-$r*$sign;
			$p=$p-$q*$sign;
			$r=$r-2*$a*$sign;
			$back=7;
			next;
		}
#step 8
		$flag=0;
		$flag++ if (($p+$q+$r+$a+$b) < -$tolerance);
		$flag++ if ((abs($p+$q+$r+$a+$b) < $tolerance) && ((2*($a+$q)+$r)>$tolerance));
		if ($flag > 0){
			$c=$a+$b+$c+$p+$q+$r;
			$p=2*$b+$p+$r;
			$q=2*$a+$q+$r;
			$back=8;
			@vec3=&sum_vv(@vec3,@vec1);
			@vec3=&sum_vv(@vec3,@vec2);
			next;
		}
		$flag=0;
	}
#determinant check
	my $determinant=&det3(@vec1,@vec2,@vec3);
	if ($determinant<0){
		@vec1=&product_vs(@vec1,-1);
		@vec2=&product_vs(@vec2,-1);
		@vec3=&product_vs(@vec3,-1);
	}
	for (my $j=0;$j<3;$j++){
		$vec1[$j]=0 if (abs($vec1[$j]) < $tolerance);
		$vec2[$j]=0 if (abs($vec2[$j]) < $tolerance);
		$vec3[$j]=0 if (abs($vec3[$j]) < $tolerance);
	}
	my @ret=(@vec1,@vec2,@vec3);
	return (@ret);
}




sub maxortho{
#cartesian ni suru
	my $b11=&product_vv(@latvec[0..2],@latvec[0..2]);
	my $b12=&product_vv(@latvec[0..2],@latvec[3..5]);
	my $b13=&product_vv(@latvec[0..2],@latvec[6..8]);
	my $b22=&product_vv(@latvec[3..5],@latvec[3..5]);
	my $b23=&product_vv(@latvec[3..5],@latvec[6..8]);
	my $y2=-($b23/$b22-$b12*$b13/$b11/$b22)/(1-$b12*$b12/$b11/$b22);
	my $y1=-($b13/$b11-$b12*$b23/$b11/$b22)/(1-$b12*$b12/$b11/$b22);
#x1,x2 kouho= x1a,x1b,y1a,y1b
	my $x1a=&ceiling($y1);
	my $x1b=&floor($y1);
	my $x2a=&ceiling($y2);
	my $x2b=&floor($y2);
	my $x1a=&ceiling($y1);
#nagasa hikaku
	my $d1=&norm(&maxortho_c($x1a,$x2a));
	my $d2=&norm(&maxortho_c($x1a,$x2b));
	my $d3=&norm(&maxortho_c($x1b,$x2a));
	my $d4=&norm(&maxortho_c($x1b,$x2b));
	my $mind=&min($d1,$d2,$d3,$d4);
	if ($d1==$mind){
		@latvec[6..8]=&maxortho_c($x1a,$x2a);
	} elsif ($d2==$mind){
		@latvec[6..8]=&maxortho_c($x1a,$x2b);
	} elsif ($d3==$mind){
		@latvec[6..8]=&maxortho_c($x1b,$x2a);
	} else{
		@latvec[6..8]=&maxortho_c($x1b,$x2b);
	}
}

sub maxortho_c{
#maxortho de c jiku no nagasa wo motomeru
#shift x1, x2 ga hikisuu, teigi wa
#http://www.csie.nuk.edu.tw/~cychen/Lattices/A%203-Dimensional%20Lattice%20Reduction%20Algorithm.pdf


	my @v1=&product_vs(@latvec[0..2],$_[0]);
	my @v2=&product_vs(@latvec[3..5],$_[1]);
	my @v3=&sum_vv(@latvec[6..8],@v1);
	my @v3=&sum_vv(@v3,@v2);
	return (@v3);
}


sub which_atom_number{
#@_ ga dono gensi to onaji ka wo shiraberu
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	unshift (@x,$_[0]);
	unshift (@y,$_[1]);
	unshift (@z,$_[2]);
	if ($dc eq "C"){;
		$num_atoms++;
		&c2d; 
		$num_atoms++;
	}
	my (@x1, @y1, @z1);
	for (my $i=0; $i<=$num_atoms; $i++){
		$x[$i]=&precise($x[$i]-&floor($x[$i]),1,1E-4);
		$y[$i]=&precise($y[$i]-&floor($y[$i]),1,1E-4);
		$z[$i]=&precise($z[$i]-&floor($z[$i]),1,1E-4);
		$x[$i]=&precise($x[$i],0,1E-4);
		$y[$i]=&precise($y[$i],0,1E-4);
		$z[$i]=&precise($z[$i],0,1E-4);
		($x1[$i],$y1[$i],$z1[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	@x=@x_orig;
	@y=@y_orig;
	@z=@z_orig;
	for (my $i=1; $i<=$num_atoms; $i++){
		my $norm=&norm(&diff_vv($x1[0],$y1[0],$z1[0],$x1[$i],$y1[$i],$z1[$i]));
		return ($i) if ($norm < 1E-3);
	}
	die ("which_atom_number: @_[0..2] to onaji ichi no gensi ga nai\n");
}

#symmetrize
sub symmetrize{
#c2d
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
#dist square
	my $dist=$_[12]*$_[12];
#image wo tsukuru
#original
	my @x1=@x;
	my @y1=@y;
	my @z1=@z;
#isometry
	my (@x2, @y2, @z2);
	for (my $i=0; $i<$num_atoms;$i++){
		my @rotate=&product_mv(@_[0..8],$x[$i],$y[$i],$z[$i]);
		($x2[$i],$y2[$i],$z2[$i])=&sum_vv(@rotate,@_[9..11]);
	}
#average
	my (@x3, @y3, @z3);
	for (my $i=0; $i<$num_atoms;$i++){
#x,y,z wo x2[j]-x1[i], etc. ni suru
#X2_rounded
		my (@x2r, @y2r, @z2r);
		for (my $j=0; $j<$num_atoms;$j++){
			my @diff;
			my @diff_round;
			@diff=&diff_vv($x2[$j],$y2[$j],$z2[$j],$x1[$i],$y1[$i],$z1[$i]);
			@diff_round=&round_array(@diff);
			($x[$j],$y[$j],$z[$j])=&diff_vv(@diff,@diff_round);
			($x2r[$j],$y2r[$j],$z2r[$j])=&diff_vv($x2[$j],$y2[$j],$z2[$j],@diff_round);
		}
#d2c de cartesian ni suru
		&d2c;
		my $min_dist=$dist;
#candidate_id
		my $min_dist_atom;
		for (my $j=0; $j<$num_atoms;$j++){
			my $norm=&product_vv($x[$j],$y[$j],$z[$j],$x[$j],$y[$j],$z[$j]);
			if ($norm < $min_dist){
				$min_dist=$norm;
				$min_dist_atom=$j;
			}
		}
#gensi ga mitsukaranai
		die ("Jyuubun chikai gensi ga mitukaranai \n") if ($min_dist == $dist);
#save average
		$x3[$i]=($x1[$i]+$x2r[$min_dist_atom])*0.5;
		$y3[$i]=($y1[$i]+$y2r[$min_dist_atom])*0.5;
		$z3[$i]=($z1[$i]+$z2r[$min_dist_atom])*0.5;
	}
#final
	@x=@x3;
	@y=@y3;
	@z=@z3;
}

#@_[0..2] ni chikai gensi bangou wo kaesu
sub near_atom{
	my $tolerance=0.01;
	for (my $i=0; $i<$num_atoms;$i++){
		my @dist=($_[0]-$x[$i], $_[1]-$y[$i], $_[2]-$z[$i]);
		@dist=&diff_vv(@dist, &round_array(@dist));
		return ($i+1) if (&norm(@dist) < $tolerance);
	}
	die ("near_atom: chikai gensi ga sonzai sinai");
}


############
#
# Phonopy
#
############


sub get_bposcar{
=pod
Phonopy:BPOSCAR wo nyusyu, yomikomu
Kuukangun jyouhou (bangou, kigou, centering, Bravais kousi,
kakutyo Bravais kousi), tengun, phonopy version type wo kaesu
@x @y @z wo uwagaki

$_[0]="noclear": phonopy no file wo nokosu

Extended Bravais symbol: Bravais lattice + type (1-3) + inversion (Y/N)
Bravais lattice and type:
cP1: #195-206
cP2: #207-230
cF1: #195-206
cF2: #207-230
cI1
tP1
tI1:c<a
tI2:c>a
oP1
oF1:a^-2>b^-2+c^-2
oF2:c^-2>a^-2+b^-2
oF3:a^-2, b^-2, c^-2 edges of triangle
oI1:c largest
oI2:a largest
oI3:b largest
oC1:a<b
oC2:a>b
oA1:b<c
oA2:b>c
hP1:#143-149, 151, 153, 157, 159-163
hP2:Other
hR1:sqrt(3)a<sqrt(2)c
hR2:sqrt(3)a>sqrt(2)c
mP1
mC1:b<asin(beta)
mC2:b>asin(beta) BZ 12-face
mC3:b>asin(beta) BZ 14-face
aP1
=cut

#kyousei D
	if ($dc eq "C"){
		$dc="D";
		&c2d;
	}
#phonopy
	&writePOSCAR ("tsubo_temp_poscar");
	&run_phonopy;

#BPOSCAR wo yomu:
#Note: motono POSCAR no cartesian genshi ichi wa keisyou sarenai!
	$scale=1;
	@latvec=@x=@y=@z=@w=();
	@latvec=&split3(`head -3 BPOSCAR | tail -1`);
	@latvec=(@latvec,&split3(`head -4 BPOSCAR | tail -1`));
	@latvec=(@latvec,&split3(`head -5 BPOSCAR | tail -1`));
	my @abcreal=&get_abc(@latvec);
#BPOSCAR ha vasp 4 to 5 ryouhou sonzai
	my $vasp5check=`head -8 BPOSCAR | tail -n 1 | cut -b-1 `;
	chop $vasp5check;
	if ($vasp5check eq "D"){
		@numspecies=&splitall(`head -7 BPOSCAR | tail -1`);
		system ("tail -n +9 <BPOSCAR>tsubo_temp_poscar");
	} else {
		@numspecies=&splitall(`head -6 BPOSCAR | tail -1`);
		system ("tail -n +8 <BPOSCAR>tsubo_temp_poscar");
	}
	$num_atoms=&numatoms;
	open BPOSCAR, "tsubo_temp_poscar";
	for (my $i=0; $i<$num_atoms; $i++){
		my $t=<BPOSCAR>;
		my @b=&split4($t);
		$x[$i]=$b[0];
		$y[$i]=$b[1];
		$z[$i]=$b[2];
	}
	close BPOSCAR;
#check phonopy version
	my @line1=&splitall(`head -1 tsubo_temp_phonopyout | tail -1`);
	$line1[1] =~ s/\'//g;
	my @phonover = split(/\./, $line1[1] );
#phonopy version type
#1.11.12 ikou ha "b"
#other "a"
	my $phonotype="b";
	if (($phonover[0] <=1)&&($phonover[1] <=10)){
		$phonotype="a" ;
	} elsif (($phonover[0] ==1)&&($phonover[1] ==11)&&($phonover[2] <12)){
		$phonotype="a" ;
	}
#kuukan gun wo check
	my ($center, $spgsymbol, $spgnumber, $ptg)=(0,0,0,0);
	if ($phonotype eq "a"){
		my @b=&splitall(`head -2 tsubo_temp_phonopyout | tail -1`);
#centering
		$center=substr($b[1], 0, 1);
#kuukan gun
		$spgsymbol=$b[1];
		$spgnumber=substr($b[2], 1, -1);
#ten gun
		@b=&split3(`head -3 tsubo_temp_phonopyout | tail -1`);
		$ptg=$b[1];
	} elsif ($phonotype eq "b"){
#centering
		my @b=&splitall(`head -2 tsubo_temp_phonopyout | tail -1`);
		$b[1] =~ s/\'//g;
		$spgsymbol=$b[1];
		$center=substr($spgsymbol, 0, 1);
#kuukan gun bangou
		@b=&splitall(`head -3 tsubo_temp_phonopyout | tail -1`);
		$spgnumber=$b[1];
#ten gun
		my @b=&splitall(`head -4 tsubo_temp_phonopyout | tail -1`);
		$b[1] =~ s/\'//g;
		$ptg=$b[1];
	} 
#inversion
	my $inversion="N";
	$inversion="Y" if ($spgnumber==2);
	$inversion="Y" if (&within($spgnumber,10,15));
	$inversion="Y" if (&within($spgnumber,47,74));
	$inversion="Y" if (&within($spgnumber,83,88));
	$inversion="Y" if (&within($spgnumber,123,142));
	$inversion="Y" if (&within($spgnumber,147,148));
	$inversion="Y" if (&within($spgnumber,162,167));
	$inversion="Y" if (&within($spgnumber,175,176));
	$inversion="Y" if (&within($spgnumber,191,194));
	$inversion="Y" if (&within($spgnumber,200,206));
	$inversion="Y" if (&within($spgnumber,221,230));
#lattice type
	my $lattice_type;
	my $bravais_symbol;
	my $bravais_symbol2="X";
	my $type=1;
	if ($spgnumber <= 2){
		$lattice_type="Triclinic";
		$bravais_symbol="aP";
	} elsif ($spgnumber <= 15){
		$lattice_type="Monoclinic";
		if ($center eq "P") {
			$bravais_symbol="mP";
		} else {
			$bravais_symbol="mS";
			$bravais_symbol2="mC";
			my $t1=$abcreal[0]*&sind($abcreal[4])/$abcreal[1];
			if ($t1 < 1){
				my $t2=1+$abcreal[0]*&cosd($abcreal[4])/$abcreal[2]-$t1*$t1;
				$type=2;
				$type=3 if ($t2 < 0);
			}
		}
	} elsif ($spgnumber <= 74){
		$lattice_type="Orthorhombic";
		if ($center eq "P"){
			$bravais_symbol="oP";
		} elsif ($center eq "I") {
			$bravais_symbol="oI";
			$type=2 if (($abcreal[0]>$abcreal[1]) && ($abcreal[0]>$abcreal[2]));
			$type=3 if (($abcreal[1]>$abcreal[0]) && ($abcreal[1]>$abcreal[2]));
		} elsif ($center eq "F") {
			$bravais_symbol="oF";
			my $a2=1/$abcreal[0]/$abcreal[0];
			my $b2=1/$abcreal[1]/$abcreal[1];
			my $c2=1/$abcreal[2]/$abcreal[2];
			if ($a2>($b2+$c2)){
				$type=1;
			} elsif ($c2>($a2+$b2)){
				$type=2;
			} else {
				$type=3;
			}
		} elsif (&within($spgnumber,38,41)) {
			$bravais_symbol="oS";
			$bravais_symbol2="oA";
			$type=2 if ($abcreal[1]>$abcreal[2]);
		} else {
			$bravais_symbol="oS";
			$bravais_symbol2="oC";
			$type=2 if ($abcreal[0]>$abcreal[1]);
		}
	} elsif ($spgnumber <= 142){
		$lattice_type="Tetragonal";
		if ($center eq "P") {
			$bravais_symbol="tP";
		} else {
			$bravais_symbol="tI";
			$type=2 if ($abcreal[2]>$abcreal[0]);
		}
	} elsif ($spgnumber <= 194){
		$lattice_type="Hexagonal_Rhombohedral";
		if ($center eq "P") {
			$bravais_symbol="hP";
			$type=2;
			$type=1 if (&within($spgnumber,143,149));
			$type=1 if ($spgnumber==151);
			$type=1 if ($spgnumber==153);
			$type=1 if ($spgnumber==157);
			$type=1 if (&within($spgnumber,159,163));
		} else {
			$bravais_symbol="hR";
			my $t=sqrt(3)*$abcreal[0]/sqrt(2)/$abcreal[2];
			$type=2 if ($t>1);
		}
	} else {
		$lattice_type="Cubic";
		if ($center eq "P"){
			$bravais_symbol="cP";
			$type=2 if (&within($spgnumber,207,230));
		} elsif ($center eq "I") {
			$bravais_symbol="cI";
		} else {
			$bravais_symbol="cF";
			$type=2 if (&within($spgnumber,207,230));
		}
	}
	$bravais_symbol2=$bravais_symbol if ($bravais_symbol2 eq "X");
	my $extended_bravais_symbol=$bravais_symbol2.$type.$inversion;
	&clear_phonopy if ($_[0] ne "noclear");
	return ($spgnumber, $spgsymbol, $center, $lattice_type, $bravais_symbol, $extended_bravais_symbol, $ptg, $phonotype);
}


sub run_phonopy{
#Phonopy wo jikkou
#Nyuuryoku file wa hikisuu, default tsubo_temp_poscar
#Syuturyoku file wa BPOSCAR. tsubo_temp_phonopyout, (PPOSCAR)
	system ("mv BPOSCAR BPOSCAR_before_tsubo ") if (-e "BPOSCAR");
	system ("mv PPOSCAR PPOSCAR_before_tsubo ") if (-e "PPOSCAR");
	my $inputfile=$_[0];
	$inputfile="tsubo_temp_poscar" if ($inputfile eq "");
	system ("phonopy --cell=tsubo_temp_poscar --symmetry --tolerance=${phonopy_tolerance} > tsubo_temp_phonopyout");
	my $t=`wc -l tsubo_temp_phonopyout`*1;
	if ($t == 2){
		&clear_phonopy;
		die ("Phonopy; Determinant of the lattice vector matrix has to be positive. \n")
	}
	while ($t == 0){
		$phonopy_tolerance*=0.5;
		if ($phonopy_tolerance < 1E-6){
			&clear_phonopy;
			die ("Phonopy; too small tolerance but segmentation fault exists. \n")
		}
		system ("phonopy --cell=tsubo_temp_poscar --symmetry --tolerance=${phonopy_tolerance} > tsubo_temp_phonopyout");
		$t=`wc -l tsubo_temp_phonopyout`*1;
	}
}

sub clear_phonopy{
#phonopy kankei no file wo sakujyo
	system ("rm tsubo_temp_poscar") if (-e "tsubo_temp_poscar");
	system ("rm tsubo_temp_phonopyout") if (-e "tsubo_temp_phonopyout");
	system ("rm BPOSCAR") if (-e "BPOSCAR");
	system ("rm PPOSCAR") if (-e "PPOSCAR");
	system ("mv BPOSCAR_before_tsubo BPOSCAR ") if (-e "BPOSCAR_before_tsubo");
	system ("mv PPOSCAR_before_tsubo PPOSCAR ") if (-e "PPOSCAR_before_tsubo");
}


#2D kuukan gun
#POSCAR no z=$_[0] no 2D kuukan gun wo kaesu
#$_[1]="renormalize": a,b kousi vector wo teisuubai site
#2D conventional cell no kousi teisuu ga 3D kousi teisuu to onaji you ni suru
#centering wa musi sareru node square to hexagonal niwa chuui
sub twoDspacegroup{
#use direct
	if ($dc eq "c"){
		$dc="d";
		&c2d;
	}
#taihi
	my @latvec_orig=@latvec;
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @w_orig=@w;
	my @namespecies_orig=@namespecies;
	my @numspecies_orig=@numspecies;
	my $num_atoms_orig=$num_atoms;
#kesu genso wo kimeru
	my @remove_list;
	my $tolerance=1e-7;
	for (my $i=$num_atoms-1; $i>=0;$i--){
		if (abs($z[$i]-$_[0]) > $tolerance){
			push (@remove_list,$i+1);
		} else {
			$z[$i]=0;
		}
	}
	&remove_atoms(@remove_list);
	die ("twoDspacegroup: z= $_[0] ni gensi ga nai! \n") if ($num_atoms == 0);

#latvec c wo kaeru

	@latvec[6..8]=&cross_vv(@latvec[0..5]);
	@latvec[6..8]=product_vs(@latvec[6..8],1E-4);
#renormalize taisaku
	my @latvec_orig2=@latvec;
	my @x_orig2=@x;
	my @y_orig2=@y;
	my @z_orig2=@z;
	my @w_orig2=@w;
	my @namespecies_orig2=@namespecies;
	my @numspecies_orig2=@numspecies;
	my $num_atoms_orig2=$num_atoms;
	my @a=&get_bposcar("noclear");
	my @xycount=(0,0);
	if ($_[1] eq "renormalize"){
		my @isometry_translation=&get_isometry_translation;
		#supercell size check
		for (my $i=0;$i<$#isometry_translation;$i+=3){
			$xycount[0]++ if (&norm(&diff_vv(@isometry_translation[$i+1],@isometry_translation[$i+2],0,0))==0);
			$xycount[1]++ if (&norm(&diff_vv(@isometry_translation[$i+0],@isometry_translation[$i+2],0,0))==0);
		}
		if (&norm(&diff_vv(@xycount,1,1)) != 0){
#rerun phonopy
			@latvec=@latvec_orig2;
			@x=@x_orig2;
			@y=@y_orig2;
			@z=@z_orig2;
			@w=@w_orig2;
			@namespecies=@namespecies_orig2;
			@numspecies=@numspecies_orig2;
			$num_atoms=$num_atoms_orig2;
			@latvec[0..2]=product_vs(@latvec[0..2],$xycount[0]);
			@latvec[3..5]=product_vs(@latvec[3..5],$xycount[1]);
			&clear_phonopy;
			@a=&get_bposcar("noclear");
		}
	}
	my @latpar=&get_abc(@latvec);
	&clear_phonopy;

	my ($sgnumber2d, $sg2d);
	if ($a[0] == 6){
		$sgnumber2d=1;
		$sg2d="p1";
	} elsif ($a[0] == 10){
		$sgnumber2d=2;
		$sg2d="p2";
	} elsif ($a[0] == 25){
		$sgnumber2d=3;
		$sg2d="pm";
	} elsif ($a[0] == 26){
		$sgnumber2d=4;
		$sg2d="pg";
	} elsif ($a[0] == 38){
		$sgnumber2d=5;
		$sg2d="cm";
	} elsif ($a[0] == 47){
		$sgnumber2d=6;
		$sg2d="p2mm";
	} elsif ($a[0] == 51){
		$sgnumber2d=7;
		$sg2d="p2mg";
	} elsif ($a[0] == 55){
		$sgnumber2d=8;
		$sg2d="p2gg";
	} elsif ($a[0] == 65){
		$sgnumber2d=9;
		$sg2d="c2mm";
#check 2 cases
#a_prim=b_prim?	
	} elsif ($a[0] == 83){
		$sgnumber2d=10;
		$sg2d="p4";
	} elsif ($a[0] == 123){
		$sgnumber2d=11;
		$sg2d="p4mm";
	} elsif ($a[0] == 127){
		$sgnumber2d=12;
		$sg2d="p4gm";
	} elsif ($a[0] == 174){
		$sgnumber2d=13;
		$sg2d="p3";
	} elsif ($a[0] == 187){
		$sgnumber2d=14;
		$sg2d="p3m1";
	} elsif ($a[0] == 189){
		$sgnumber2d=15;
		$sg2d="p31m";
	} elsif ($a[0] == 175){
		$sgnumber2d=16;
		$sg2d="p6";
	} elsif ($a[0] == 191){
		$sgnumber2d=17;
		$sg2d="p6mm";
	} else {
		die ("twoDspacegroup subroutine de 3D kuukangun bangou ga okasii: 3D kuukangun @a[0..1] \n");
	}

#modosu
	@latvec=@latvec_orig;
	@x=@x_orig;
	@y=@y_orig;
	@z=@z_orig;
	@w=@w_orig;
	@namespecies=@namespecies_orig;
	@numspecies=@numspecies_orig;
	$num_atoms=$num_atoms_orig;
	return ($sgnumber2d,$sg2d,@latpar[0..2],@xycount);
}

#Phonopy:BPOSCAR wo moto ni kessyougaku/standard 
#conventional/primitive cellwo motomeru
#saisyo ni bposcar wo tsukuru (@latvec wo crystallographic conventional ni suru
#$_[0] ha "conventional" or "primitive"
#Direct de gensi ichi wo ataeru, direct de kaesu

#Kessyougaku version
sub get_crystallographic_cell{
	my @t=&get_bposcar;
	my $center=$t[2];
	my $lattice_type=$t[3];
#Extended Bravais no saisyo no ni moji
	my $bravais=substr($t[5],0,2);
#lattice parameters real
	my @abcreal=&get_abc(@latvec);
#new primitive lattice vector
	my @latvec_new;
#default primitive
	my $default="P";
#jiku wo torinaosu, cartesian ni suru
	&d2c; 
#primitive
	if ($_[0] eq "primitive"){
		if ($center eq "P"){
			&get_crystallographic_cell_get_latvec("P");
		} elsif ($center eq "I"){
			&get_crystallographic_cell_get_latvec("I");
		} elsif ($center eq "F"){
			&get_crystallographic_cell_get_latvec("F");
		} elsif ($center eq "R"){
			&get_crystallographic_cell_get_latvec("R");
		} elsif ($bravais eq "oC"){
			&get_crystallographic_cell_get_latvec("oC");
		} elsif ($bravais eq "oA"){
			&get_crystallographic_cell_get_latvec("oA");
		} elsif ($bravais eq "mC"){
			&get_crystallographic_cell_get_latvec("mC");
		} else {
			die ("bug in get_crystallographic_cell\n");
		}
	} else {
		&get_crystallographic_cell_get_latvec("P");
	}
#return extended Bravais
	return ($t[5]);
}


# kousi teisuu wo kaku
# P, C, I, F, R no latvec wo kaesu
# centering wo hikisuu to site ataeru
# @latvec wo tsukau
# atarasii latvec wo kaesu
# zahyo ha cartesian to suru! direct de kaesu
sub get_crystallographic_cell_get_latvec{
#standardize
	&c2d;
	my @abcreal=&get_abc(@latvec);
	my @t=&clat(@abcreal,"s");
	@latvec=($t[0],0,0,$t[1],$t[2],0,@t[3..5]);
	if ($_[0] ne "P"){
		&d2c;
		if ($_[0] eq "F"){
			@latvec=&product_mm(0,1,1,1,0,1,1,1,0,@latvec);
			@latvec=&product_vs(@latvec,1/2);
		} elsif ($_[0] eq "I"){
			@latvec=&product_mm(-1,1,1,1,-1,1,1,1,-1,@latvec);
			@latvec=&product_vs(@latvec,1/2);
		} elsif ($_[0] eq "R"){
			@latvec=&product_mm(2,1,1,-1,1,1,-1,-2,1,@latvec);
			@latvec=&product_vs(@latvec,1/3);
		} elsif ($_[0] eq "oC"){
			@latvec=&product_mm(1,-1,0,1,1,0,0,0,2,@latvec);
			@latvec=&product_vs(@latvec,1/2);
		} elsif ($_[0] eq "oA"){
			@latvec=&product_mm(0,1,-1,0,1,1,2,0,0,@latvec);
			@latvec=&product_vs(@latvec,1/2);
		} elsif ($_[0] eq "mC"){
			@latvec=&product_mm(1,1,0,-1,1,0,0,0,2,@latvec);
			@latvec=&product_vs(@latvec,1/2);
		} else {
			die ("bug in get_crystallographic_cell_get_latvec\n");
		}
		&c2d;
		&unique();
	}
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
}


#Hinuma reduced (triclinic only)
#Direct de zahyo wo uketoru, Direct de kaesu
#Acute ka Obtuse mo kaesu
# get niggli-reci poscar
sub get_hinuma_reduced{
	my @latvec_reci=&tmatrix3(&invmatrix3(@latvec));
	my @latvec_reci_niggli=&get_niggli_latvec(@latvec_reci);
	my @latvec_niggli=&tmatrix3(&invmatrix3(@latvec_reci_niggli));
	my @conversion_matrix=&round_array(&product_mm(@latvec,&invmatrix3(@latvec_niggli)));
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i], @conversion_matrix);		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	@latvec=@latvec_niggli;
#Direct to cartesian
	&d2c;
#ka*kb*k_gamma wo 90 do ni chikazukeru
	my @abcreal=&get_abc(@latvec);
	my @recilatvec=&reci_latvec(@latvec);
	my @abcreci=&get_abc(@recilatvec);
	my $diffx=abs($abcreci[1]*$abcreci[2]*&cosd($abcreci[3]));
	my $diffy=abs($abcreci[2]*$abcreci[0]*&cosd($abcreci[4]));
	my $diffz=abs($abcreci[0]*$abcreci[1]*&cosd($abcreci[5]));
	if ($diffx == &min($diffx,$diffy,$diffz)){
		&rotate_axis("bca");
	} elsif ($diffy == &min($diffx,$diffy,$diffz)){
		&rotate_axis("cab");
	} 
	@abcreal=&get_abc(@latvec);
	@recilatvec=&reci_latvec(@latvec);
	@abcreci=&get_abc(@recilatvec);
#zenbu >90 matawa <90 ni suru
	if (($abcreci[3] <= 90) && ($abcreci[4] > 90) && ($abcreci[5] > 90)){
		&rotate_axis("a-b-c");
	} elsif (($abcreci[3] > 90) && ($abcreci[4] <= 90) && ($abcreci[5] > 90)){
		&rotate_axis("-ab-c");
	} elsif (($abcreci[3] > 90) && ($abcreci[4] > 90) && ($abcreci[5] <= 90)){
		&rotate_axis("-a-bc");
	} elsif (($abcreci[3] > 90) && ($abcreci[4] <= 90) && ($abcreci[5] <= 90)){
		&rotate_axis("a-b-c");
	} elsif (($abcreci[3] <= 90) && ($abcreci[4] > 90) && ($abcreci[5] <= 90)){
		&rotate_axis("-ab-c");
	} elsif (($abcreci[3] <= 90) && ($abcreci[4] <= 90) && ($abcreci[5] > 90)){
		&rotate_axis("-a-bc");
	}
	@abcreal=&get_abc(@latvec);
	@recilatvec=&reci_latvec(@latvec);
	@abcreci=&get_abc(@recilatvec);
	my $type="Obtuse";
	$type="Acute" if ($abcreci[5] < 90);
	&get_crystallographic_cell_get_latvec("P");
	return $type;
}


sub get_isometry{
#Phonopy:Isometry no list wo tsukuru: rotationx9 + translationx3 wo kurikaesi
#Genzai no POSCAR data wo siyou
#space group number wo saisyo ni kaesu
	&writePOSCAR ("tsubo_temp_poscar");
	my @run_phonopy_data=&run_phonopy if ((! -e "tsubo_temp_phonopyout") || (! -e "tsubo_temp_poscar"));
	my @isometry;
#first rotation line
	my $firstrot=`grep -n "rotation:" tsubo_temp_phonopyout | head -1 | cut -d ":" -f 1`*1;
	system ('tail -n +'.${firstrot}.' <tsubo_temp_phonopyout | sed -e "s/\[\|\]\|,//g" > tsubo_temp_poscar');

	open TEMP, "tsubo_temp_poscar";
	my @a=&splitall(my $t=<TEMP>);
	while ($a[1] eq "rotation:"){
		for (my $i=0;$i<4;$i++){
			@a=&splitall(my $t=<TEMP>);
			push(@isometry,@a[1..3]);
		}
		@a=&splitall(my $t=<TEMP>);
	}
	close TEMP;

#kuukan gun bangou: 
#	my @b=&split3(`head -2 tsubo_temp_phonopyout | tail -1`);
	&clear_phonopy;
	return (@isometry);
}



#Real Isometry 
#tadasii isometry wo kaesu (vector part no precision wo ageru)
#nyuuryoku ha (isometry)
sub isometry_real{
#tadasii isometry wo nyuusyu
#site (gensi 1)
	my @site=($x[0],$y[0],$z[0]);
	my @Wx=&product_mv(@_[0..8],@site);
	my @tempsite=&sum_vv(@Wx,@_[9..11]);
	my $imageatom=&which_atom_number(@tempsite);
#tadasii vector part
	my @a=&diff_vv($x[$imageatom-1],$y[$imageatom-1],$z[$imageatom-1],@Wx);
	return (@_[0..8],@a);
}


sub get_isometry_translation{
#Phonopy:Translation no list wo tsukuru: translationx3 wo kurikaesi
#Matrix ga identity nomi tyuusyutu
#Genzai no POSCAR data wo siyou
	my @all_isometry=&get_isometry;
#matrix ga identity no vector bubun wo nyuusyu
	my @translation;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		if (&norm(&diff_vv(@all_isometry[$i..($i+8)],1,0,0,0,1,0,0,0,1)) == 0){
			push (@translation, @all_isometry[($i+9)..($i+11)]);
		}
	}
	return (@translation);
}


sub get_surface_isometry{

#Non polar surface ni tsukaeru isometry wo sagasu
#Gutaiteki niha (*,*,*,*,*,*,0,0,-1) no katachi
#Gaitou isometry wo hairetsu to shite kaesu 
	my @list;
	my @all_isometry=&get_isometry;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		my @a=@all_isometry[$i..($i+11)];
#tsukaeruka hantel
		push (@list, @a) if (&norm(&diff_vv(@a[6..8],0,0,-1)) == 0);
	}
	return @list;
}


#Isometry no shurui wo kaesu
#type (-1, 4, m nado), glide/screw/translation part/center(inversion/rotoinversion), plane/axis info
#axis info: a..l ga syuturyoku
# {{a,b,c},{d,e,f},{g,h,i}}*{x,y,z}T={j,k,l}T ni naru x,y,z ga axis/plane
#hutei ga detekuru renritu houteisiki dakara toku noga muzukasii
#proper no baai ha glide/screw part ha 0 0 0 or lattice vector
#supercell no baai ha translation vector ga jimei de nai node chuui
#Isometry (suuji x 12) ga nyuuryoku
sub get_isometry_type{
	my @return;
	my $flag;
	my $trace=$_[0]+$_[4]+$_[8];
	my $det=&det3(@_[0..8]);
#type
	if (($trace==-1)&&($det==1)){
		$return[0]=2;
		$flag=2;
	} elsif (($trace==0)&&($det==1)){
		$return[0]=3;
		$flag=3;
	} elsif (($trace==1)&&($det==1)){
		$return[0]=4;
		$flag=4;
	} elsif (($trace==2)&&($det==1)){
		$return[0]=6;
		$flag=6;
	} elsif (($trace==3)&&($det==1)){
		$return[0]=1;
		$flag=1;
	} elsif (($trace==-3)&&($det==-1)){
		$return[0]=-1;
		$flag=-1;
	} elsif (($trace==-2)&&($det==-1)){
		$return[0]=-6;
		$flag="rotoinversion";
	} elsif (($trace==-1)&&($det==-1)){
		$return[0]=-4;
		$flag="rotoinversion";
	} elsif (($trace==0)&&($det==-1)){
		$return[0]=-3;
		$flag="rotoinversion";
	} elsif (($trace==1)&&($det==-1)){
		$return[0]="m";
		$flag="mirror";
	} else {
		die ("get_isometry_type: determinant to trace no kumiawase ga okasii: @_[0..8] trace $trace determinant $det\n");
	}
#symmetry operation: 1
	if ($flag eq "1"){
		@return[1..15]=($_[9]*0.5,$_[10]*0.5,$_[11]*0.5,"na","na","na","na","na","na","na","na","na","na","na","na");
#symmetry operation: -1
	} elsif ($flag eq "-1"){
		@return[1..15]=($_[9]*0.5,$_[10]*0.5,$_[11]*0.5,"na","na","na","na","na","na","na","na","na","na","na","na");
#symmetry operation: rotoinversion
	} elsif ($flag eq "rotoinversion"){
#get center: (W,w)=x
		@return[1..3]=&product_mv(&invmatrix3(&sum_vv(@_[0..8],-1,0,0,0,-1,0,0,0,-1)),&product_vs(@_[9..11],-1));
#get axis: symmetry operation (W,w)^2
		my @symop=&product_symop(@_[0..11],@_[0..11]);
#check matrix, hutei ga 1 aru hazu
#dekireba 1 ni oshituketai ga y ka z ga koritu siteitara muri
		@return[4..15]=(&sum_vv(@symop[0..8],-1,0,0,0,-1,0,0,0,-1),@symop[9..11]);
	} elsif ($flag eq "mirror"){
#glide: translation part
		my @symop=&product_symop(@_[0..11],@_[0..11]);
		@return[1..3]=&product_vs(@symop[9..11],0.5);
		@return[4..12]=&sum_vv(@_[0..8],-1,0,0,0,-1,0,0,0,-1);				@return[13..15]=&diff_vv(@_[9..11],@return[1..3]);
#rotation
	} else {
		my @symop=@_[0..11];
		for (my $i=1; $i<$flag; $i++){
			@symop=&product_symop(@symop,@_[0..11]);
		}
		@return[1..3]=&product_vs(@symop[9..11],1/$flag);
		@return[4..12]=(&sum_vv(@_[0..8],-1,0,0,0,-1,0,0,0,-1));
		@return[13..15]=&diff_vv(@_[9..11],@return[1..3]);
	}

	return (@return);
}

#Isometry image 
#isometry ni taisi gensi ka zahyou no image wo kaesu
#nyuuryoku ha (isometry, gensi bangou ) OR (isometry, x,y,z)
sub isometry_image{
#flag 0: genshi
	my $flag=0;
#flag 1: xyz
	$flag=1 if ($_[13] ne "");
#tadasii isometry wo nyuusyu
#site
	my @site;
	if ($flag == 0){
		@site=($x[$_[12]-1],$y[$_[12]-1],$z[$_[12]-1]);
	} else {
		@site=($x[0],$y[0],$z[0]);
	}
	my @Wx=&product_mv(@_[0..8],@site);
	my @tempsite=&sum_vv(@Wx,@_[9..11]);
	my $imageatom=&which_atom_number(@tempsite);
	return ($imageatom) if ($flag == 0);
#tadasii vector part
	my @realw=&diff_vv($x[$imageatom-1],$y[$imageatom-1],$z[$imageatom-1],@Wx);
	my @image=&sum_vv(&product_mv(@_[0..8],@_[12..14]),@realw);
	for (my $i=0;$i<3;$i++){
		$image[$i]=&precise($image[$i]-&floor($image[$i]),1,1E-9);
		$image[$i]=&precise($image[$i],0,1E-9);
	}
	my @image=&site_in_cell(@image);
	return @image;
}


#lattice point ga genshi ichi to naru empty cell wo tsukuru
# hikisuu $h,$k,$l OR @transformation_matrix
# @x @y @z @w @numspecies @namespecies wo hakai suru node chuui!!
# hennkan ni tsukatta gyouretu wo kaesu
sub make_empty_cell{
#translation wo kakuho
	my @transformation_matrix=@_[0..8];
	if ($_[3] eq ""){
#hkl-mode
		@transformation_matrix=&hkl_transformation_matrix(@_[0..2]);
	} else {
#T-matrix = @_ mode
		die ("make_empty_cell: transformation_matrix no bubun no hikisuuga 9 ko nai : @transformation_matrix\n") if ($_[8] eq "");
		die ("make_empty_cell: transformation_matrix no determinant <= 0") if (&det3(@transformation_matrix)<=0);
	}

#empty cell no jyunbi
	my @translation=&get_isometry_translation;
	@x=@y=@z=@w=();
	@numspecies=(0);
	@namespecies=("X");
#empty cell wo tsukuru
#translation wa D nano de zahyou ha subete D
#koko de @x @y @z @w ga kimaru
	for (my $i=0; $i<$#translation; $i+=3){
		$numspecies[0]++;
		push (@x,$translation[$i]);
		push (@y,$translation[$i+1]);
		push (@z,$translation[$i+2]);
		push (@w,"X");
	}
#sanitize coordinate
	$num_atoms=&numatoms;
	for (my $i=0; $i<$num_atoms; $i++){
		$x[$i]=&round($x[$i]*$num_atoms)/$num_atoms;
		$y[$i]=&round($y[$i]*$num_atoms)/$num_atoms;
		$z[$i]=&round($z[$i]*$num_atoms)/$num_atoms;
	}
#kari no supercell wo tsukuru
	&get_supercell (@transformation_matrix);
# $num_atoms kousin
	$num_atoms=&numatoms;
}




########
#
# Surface related
#
########


#slab wo tsukuru
sub make_slab{
#input: minz,maxz,shiftz,tolerance
#tolerance
#slab boundary wo aimai ni suru
#writeposcar: "no" de poscar wo kakanai

#Direct ni suru
	&c2d if ($dc eq "C");
#shoki settei
	my ($minz, $maxz, $shiftz, $tolerance)=@_[0..3];
	$minz=0 if ($_[0] eq "");
	$maxz=1-1E-6 if ($_[1] eq "");
	$shiftz=0 if ($_[2] eq "");
	$tolerance=0 if ($_[3] eq "");
	die ("make_slab $maxz - $minz + 2 * $tolerance >= 1: quit\n") if (($maxz-$minz+2*$tolerance>=1));
	die ("make_slab $maxz <= $minz : quit\n") if ($maxz<=$minz);
#ikkain minz=0 ni z wo shift = -($minz-$tolerance)
	$num_atoms=&numatoms;
	for (my $i=0; $i<$num_atoms; $i++){
		$z[$i]-=($minz-$tolerance);
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	my (@x_new, @y_new, @z_new, @w_new);
	my @numspecies_new=&diff_vv(@numspecies,@numspecies);
#z no kakunin
	my $count=0;
	for (my $i=0;$i<=$#numspecies;$i++){
		my (@x_temp, @y_temp, @z_temp, @w_temp);
#hituyouna gensi wo tuika
		for (my $j=0; $j<$numspecies[$i]; $j++){
#maxz+tolerance-(minz-tolerance)
			if (($z[$count] >= 0) && ($z[$count] <= ($maxz-$minz+2*$tolerance))){
				push (@x_temp, $x[$count]);
				push (@y_temp, $y[$count]);
#awasete +shiftz
				push (@z_temp, $z[$count]+$shiftz);
				push (@w_temp, $w[$count]);
				$numspecies_new[$i]++;
			}
			$count++;
		}
		my @order = sort { $z_temp[$a] <=> $z_temp[$b] or $x_temp[$a] <=> $x_temp[$b] or $y_temp[$a] <=> $y_temp[$b]  } 0 .. $#z_temp;
		@x_new = (@x_new, @x_temp[@order]);
		@y_new = (@y_new, @y_temp[@order]);
		@z_new = (@z_new, @z_temp[@order]);
		@w_new = (@w_new, @w_temp[@order]);
	}
#gensisuu wo update
	@numspecies=@numspecies_new;
	$num_atoms=&numatoms;
	@x=@x_new;
	@y=@y_new;
	@z=@z_new;
	@w=@w_new;
#site wo modosu
	for (my $i=0; $i<$num_atoms; $i++){
		$z[$i]+=($minz-$tolerance);
	}
	&d2c if ($dc eq "C");
}



#HLK houkou no supercell wo nyusyu : PRIMITIVE or SUPER cell
#hikisuu wa (supercell OR primitive) H K L [asis]

#ALGORITHM MINOR CHANGE 2679/6/7 (use make_empty_cell)
sub get_hkl_cell{
	my @latvec_initial=@latvec;
	my $dcorig=$dc;
#force d
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
#conventional cell ni suru (hituyou nara)
	&get_bposcar if ($_[4] ne "asis");
#original cell wo taihi: C de hozon
	&d2c;
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @w_orig=@w;
	my @latvec_orig=@latvec;
	my @numspecies_orig=@numspecies;
	my @namespecies_orig=@namespecies;

#empty cell wo tsukuru
#empty cell no zahyou ha D (motomoto D dakedo)!
	&c2d;
	$dc="D";
	&make_empty_cell(@_[1..3]);
# hituyou na isometry wo kakutei
# In-plane w1 (*,0,0): 
	my $w11=1 ;
	for (my $i=0; $i<$num_atoms; $i++){
		if (($y[$i] == 0) && ($z[$i] == 0) && ($x[$i] > 0)){
			$w11=$x[$i] if ($x[$i] < $w11);
		}
	}
# In-plane w2 (*,*,0)	
	my $w21=0 ;
	my $w22=1 ;
	for (my $i=0; $i<$num_atoms; $i++){
		if (($z[$i] == 0) && ($y[$i] > 0)){
			if (($y[$i] < $w22) && ($x[$i] < $w11)){
				$w22=$y[$i];
				$w21=$x[$i];
			}
		}
	}
# Out-of-plane w3 (0,0,*)
	my $w31=0 ;
	my $w32=0 ;
	my $w33s=1 ;
	for (my $i=0; $i<$num_atoms; $i++){
		if (($x[$i] == 0) && ($y[$i] == 0) && ($z[$i] > 0)){
#Supercell mode
			$w33s=$z[$i] if ($z[$i] < $w33s);
		}
	}
#C ni suru
	&d2c;
#latvec wo super ni suru
	@latvec=&product_mm($w11,0,0,$w21,$w22,0,0,0,$w33s,@latvec);
	my @latvec_temp=&gaussian_reduction_2d(@latvec[0..5],0,0,0,0,0,0,0);
	@latvec[0..5]=@latvec_temp[0..5];
	@latvec[3..5]=&product_vs(@latvec[3..5],-1) if (&det3(@latvec)<0);
#super wo hozon

	my @latvec_super=@latvec;
#D ni suru -- primitive sagashi

	&c2d;
#unique + in cell
	&unique();
# Out-of-plane w3 (*,*,*)
	my $w31=0 ;
	my $w32=0 ;
	my $w33=1 ;
	for (my $i=0; $i<$num_atoms; $i++){
		if (($z[$i] < $w33) && ($z[$i] > 0)){
			$w31=$x[$i];
			$w32=$y[$i];
			$w33=$z[$i];
		}
	}
	@latvec=&product_mm(1,0,0,0,1,0,$w31,$w32,$w33,@latvec);
#moto no kousi ni modosu (c de hozon)
	@x=@x_orig;
	@y=@y_orig;
	@z=@z_orig;
	@w=@w_orig;
	@numspecies=@numspecies_orig;
	@namespecies=@namespecies_orig;	
	$num_atoms=&numatoms;
# D ni suru, site unique
	&c2d;
	&unique();

############ 2680/11 yori mae dewa ############
############ koko wo comment       ############
	&d2c;
	&maxortho;
	&c2d;
###############################################


# supercell?
	if ($_[0] eq "supercell"){
# transformation matrix
		my @matrix=&product_mm(@latvec_super,&invmatrix3(@latvec));
		@matrix=&round_array(@matrix);
		&get_supercell(@matrix);
	}
	if ($dcorig eq "C"){
		$dc="C";
		&d2c;
	}
}


sub hkl_transformation_matrix{
#hkl-mode
#$_[3]="nonCMS" no baai, h, k, l no uchi 2 ko ga 0 no matrix ga Comp Mater Sci 113 221 to kotonaru
	my ($h, $k, $l)=@_[0..2];
	my @t;
	die ("H K L no H ga seisuu de nai: $h \n") if ($h !~ /^[-]?[0-9]+$/);
	die ("H K L no K ga seisuu de nai: $k \n") if ($k !~ /^[-]?[0-9]+$/);
	die ("H K L no L ga seisuu de nai: $l \n") if ($l !~ /^[-]?[0-9]+$/);
	die ("H K L ga zenbu 0 \n") if (&norm(&diff_vv($h,$k,$l,0,0,0)) == 0);
#make transformation matrix
	if (($k==0) && ($l==0)){
		@t=(0,1,0,0,0,1,1,0,0);
		@t=(0,1,0,0,0,1,$h,0,0) if ($_[3] eq "nonCMS");
	} elsif (($h==0) && ($l==0)){
		@t=(0,0,1,1,0,0,0,1,0);
		@t=(0,0,1,1,0,0,0,$k,0) if ($_[3] eq "nonCMS");
	} elsif (($h==0) && ($k==0)){
		@t=(1,0,0,0,1,0,0,0,1);
		@t=(1,0,0,0,1,0,0,0,$l) if ($_[3] eq "nonCMS");
	} elsif ($h==0){
		@t=(1,0,0,0,$l,-$k,0,$k,$l);
	} elsif ($k==0){
		@t=(0,1,0,$l,0,-$h,$h,0,$l);
	} elsif ($l==0){
		@t=(0,0,1,$k,-$h,0,$h,$k,0);
	} else {
		@t=($k,-$h,0,$l,0,-$h,$h,$k,$l);
	}
	@t=&product_mm(1,0,0,0,-1,0,0,0,1,@t) if (&det3(@t)<0);
	return (@t);
}


sub nonpolar_area_sort{
#Nonpolar indices (h1, k1, l1, h2, l2, k2...) wo hikisuu
#Hikisuu wa 3 no baisuu, 0 igai
#Cell ha conventional wo katei
#saisyou no $_[0] (cap) bai ijyou wa hyouji sinai
	my (@h, @k, @l, @area);
	my $count=($#_)/3;
#Empty lattice wo tsukuru
	$dc="D";
#hairetu wo yomu 
	for (my $i=0;$i<$count;$i++){
		$h[$i]=$_[3*$i+1];
		$k[$i]=$_[3*$i+2];
		$l[$i]=$_[3*$i+3];
	}
#cell wo backup
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @w_orig=@w;
	my @latvec_orig=@latvec;
	my @numspecies_orig=@numspecies;
	my @namespecies_orig=@namespecies;

#menseki wo nyuusyu
	for (my $i=0;$i<$count;$i++){
#cell wo modosu
		@x=@x_orig;
		@y=@y_orig;
		@z=@z_orig;
		@w=@w_orig;
		@latvec=@latvec_orig;
		@numspecies=@numspecies_orig;
		@namespecies=@namespecies_orig;	
		$num_atoms=&numatoms;
		&get_hkl_cell("primitive", $h[$i], $k[$i], $l[$i]);
		$area[$i]=&norm(&cross_vv(@latvec[0..5]));
	}
#sort
	my @order_key=&order_key_number(@area);
	@h=@h[@order_key];
	@k=@k[@order_key];
	@l=@l[@order_key];
	@area=@area[@order_key];

	for (my $i=0;$i<$count;$i++){
		my $t=$area[$i]/$area[0];
		if ($t <= $_[0]){
			print ("$h[$i] $k[$i] $l[$i] ");
			printf ("%2.2f %2.2f\n", $area[$i], $area[$i]/$area[0]);
		}
	}
	exit;
}

sub get_surface_slab_info{
#Surface slab no center wo nyuusyu
#Syuki no gyakusuu to center wo kaesu
	my @list;
	my @surface_isometry=&get_surface_isometry;
#wz no list wo tyuushutu
	for (my $i=0;$i<$#surface_isometry;$i+=12){
		push (@list,$surface_isometry[($i+11)]);
	}
	@list=&unique_array(@list);
	return ("Polar") if ($list[0] eq "");
	my $period=($#list+1)*2;
	my $shift=&round(&min(@list)*24)/48;
	return ($period, $shift);
}



#surface slab wo tsukuru
#normal: imano directory ni ireru
#directory: h_k_l directory wo tsukutte ireru
#hikisuu: ("normal" or "directory") slab_min_thick vac_min_thick h k l [asis]
#NONPOLAR A,B,C!
sub make_surface_slab{
	my $make_slab_tolerance=1E-8;
#force d
	if ($dc eq "C"){
		$dc="D";
		&c2d;
	}
	my $dir=${_[3]}."_".${_[4]}."_".${_[5]};
#get data
	my $slab_min_thick=$_[1];
	my $vac_min_thick=$_[2];
	my @termination_data=&get_termination_polarity(@_[3..6]);

#modottekuru data ha flag, termination type, unit thickness, numcenter, shift...
	shift @termination_data;
# @termination_type: 1,2,3,4 for type A, B (tyouhuku=onaji type)
	my @termination_type=splice(@termination_data,0,4);
	system ("head -1 temp_tsubo_polarity_data ");
	system ("head -4 temp_tsubo_polarity_data | tail -2");
	system ("rm temp_tsubo_polarity_data");
	print ("Atsusa data \n");
#cell atsusa
	my $prim_thick=(shift (@termination_data))*1;
	my $unit_thick=$prim_thick*2/(shift (@termination_data));
	my $shift=shift (@termination_data);
	printf ("Saisyou slab thickness %2.2f Saisyou vacuum thickness %2.2f \n", $slab_min_thick, $vac_min_thick);
	printf ("Primitive cell thickness %2.2f unit thickness %2.2f shift zPc %0.2f \n", $prim_thick, $unit_thick,$shift);
	while ($termination_data[0] ne ""){
#filename
		my $filename;
#half-integer=1, integer=2
		my $integer=shift (@termination_data);
#center shift no=0, yes=1
		my $shiftcenter=shift (@termination_data);
#slab model type (1-4)
		my $modeltype=($integer-1)*2+$shiftcenter;
#lower boundary
		my $lower_boundary=shift (@termination_data);
#Termination
		my $termination_id=$termination_type[$modeltype];
		my $termination_type=shift(@termination_data);

		if ($termination_type eq "Type_A"){
			$termination_id="E" if ($termination_id eq "1");
			$termination_id="F" if ($termination_id eq "2");
			$termination_id="G" if ($termination_id eq "3");
			$termination_id="H" if ($termination_id eq "4");
			system ("mkdir $dir") if (($_[0] eq "directory") && (! -e $dir));
			print ("Nonpolar type A surface: ");
			$filename="POSCAR.NPA.".$modeltype.".".$termination_id.".";
			$filename=$dir."/".$filename if ($_[0] eq "directory");
		} elsif ($termination_type eq "Type_B"){
			$termination_id="J" if ($termination_id eq "1");
			$termination_id="K" if ($termination_id eq "2");
			$termination_id="L" if ($termination_id eq "3");
			$termination_id="M" if ($termination_id eq "4");
			system ("mkdir $dir") if (($_[0] eq "directory") && (! -e $dir));
			print ("Nonpolar type B surface: ");
			$filename="POSCAR.NPB.".$modeltype.".".$termination_id.".";
			$filename=$dir."/".$filename if ($_[0] eq "directory");
		} elsif ($termination_type eq "Type_C"){
			$termination_id="U" if ($termination_id eq "5");
			$termination_id="V" if ($termination_id eq "6");
			$termination_id="W" if ($termination_id eq "7");
			$termination_id="X" if ($termination_id eq "8");
			system ("mkdir $dir") if (($_[0] eq "directory") && (! -e $dir));
			print ("Nonpolar type C surface: ");
			$filename="POSCAR.NPC.".$modeltype.".".$termination_id.".";
			$filename=$dir."/".$filename if ($_[0] eq "directory");
		} else {
			die ("termination_type in make_surface_slab ga okasii: $termination_type \n");
		}
#slab thickness
		my $slab_thick;
		if ($integer==1){
#half-integer
			$slab_thick=(&ceiling($slab_min_thick/$unit_thick-0.5)+0.5)*$unit_thick;
		} else {
#integer
			$slab_thick=&ceiling($slab_min_thick/$unit_thick)*$unit_thick;
		}
#cell thickness
		my $n=&ceiling(($slab_thick+$vac_min_thick)/$prim_thick);
		my $cell_thick=$n*$prim_thick;
		printf ("Slab thickness %2.2f(%0.2f) Cell thickness %2.2f N $n ",$slab_thick,$slab_thick/$cell_thick,$cell_thick);
		my $zminus=$lower_boundary/$n;
		my $zplus=$zminus+$slab_thick/$cell_thick;
#makeslab
		printf ("z- %0.6f z+ %0.6f \n",$zminus,$zplus);
		my @x_orig=@x;
		my @y_orig=@y;
		my @z_orig=@z;
		my @w_orig=@w;
		my @latvec_orig=@latvec;
		my @numspecies_orig=@numspecies;
		my $makeslab_return=0;
#make supercell,slab
		if ($termination_type eq "Type_C"){
			print ("Termination $termination_id wo kentou \n");
			$makeslab_return=&reconstruct_typeC_surface($n,$zminus,$zplus,$_[6]);
		} else {
			&get_supercell(1,1,$n);
			&make_slab($zminus-$make_slab_tolerance,$zplus+$make_slab_tolerance);
		}
		if ($makeslab_return != 999){
#c-axis wo naosu:  Type C de slab ga tsukurenai baai wa pass
			&d2c;
			&maxortho;
			&c2d;
			for (my $i=0; $i<$num_atoms; $i++){
				($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
			}
			my $z_length=int(&norm(@latvec[6..8]));
			$filename=$filename.$num_atoms.".".$z_length;
			&writePOSCAR ($filename);
		}
#modosu
		@x=@x_orig;
		@y=@y_orig;
		@z=@z_orig;
		@w=@w_orig;
		@latvec=@latvec_orig;
		@numspecies=@numspecies_orig;
		$num_atoms=&numatoms;


#alternative @C
		if (($makeslab_return > 0) && ($makeslab_return != 999)){
			print ("make alternative \n");
#taihi
			@x_orig=@x;
			@y_orig=@y;
			@z_orig=@z;
			@w_orig=@w;
			@latvec_orig=@latvec;
			@numspecies_orig=@numspecies;
			$termination_id="Q" if ($termination_id eq "U");
			$termination_id="T" if ($termination_id eq "V");
			$termination_id="Y" if ($termination_id eq "W");
			$termination_id="Z" if ($termination_id eq "X");
			$makeslab_return=&reconstruct_typeC_surface($n,$zminus,$zplus,"Y");
			if ($makeslab_return != 999){
		#c-axis wo naosu:  Type C de slab ga tsukurenai baai wa pass
				&d2c;
				&maxortho;
				&c2d;
				for (my $i=0; $i<$num_atoms; $i++){
					($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
				}
				my $z_length=int(&norm(@latvec[6..8]));

				$filename="POSCAR.NPC.".$modeltype.".".$termination_id.".";
				$filename=$dir."/".$filename if ($_[0] eq "directory");

				$filename=$filename.$num_atoms.".".$z_length;
				&writePOSCAR ($filename);
			}
#modosu
			@x=@x_orig;
			@y=@y_orig;
			@z=@z_orig;
			@w=@w_orig;
			@latvec=@latvec_orig;
			@numspecies=@numspecies_orig;
			$num_atoms=&numatoms;
		}
	}
}

sub get_inplane_lattice_vector{
#hkl slab no inplane lattice vector wo sagasu
#vector no nagasa to kakudo wo kyoku zahyou de kaesu

#empty cell ni suru
#cell wa hakai suru
#force d
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
#conventional cell ni suru (hituyou nara)
	&get_bposcar if ($_[3] ne "asis");
#empty cell no zahyou ha D (motomoto D dakedo)!
#z=0 nomi hokan
	my @latvec_first=@latvec;
	&make_empty_cell(@_[0..2]);
	my @latvec_rotated=@latvec;
	for (my $i=$num_atoms-1; $i>=0; $i--){
		if (&precise(@z[$i],0)!=0){
			splice(@x,$i,1);
			splice(@y,$i,1);
			splice(@z,$i,1);
			$num_atoms--;
		} else { 
			$z[$i]=0;
		}
	}
	$numspecies[0]=$num_atoms;
#sousaku hankei
	my $range=3;
	my $h=&ceiling($range*&product_vv(@latvec[0..2],@latvec[0..2])/&norm(&cross_vv(@latvec[0..5])));
#supercell
	&get_supercell($range*2,$h*2,1);

#shift origin, restore latvec
	for (my $i=$num_atoms-1; $i>=0; $i--){
		$x[$i]-=0.5;
		$y[$i]-=0.5;
#genten wo korosu
		if ((&precise($x[$i],0)==0) && (&precise($y[$i],0)==0)){
			splice(@x,$i,1);
			splice(@y,$i,1);
			splice(@z,$i,1);
			$num_atoms--;
		}
	}
	&d2c;
	@latvec=@latvec_first;
	my @xc=@x;
	my @yc=@y;
	my @zc=@z;
	&c2d;
	for (my $i=0; $i<$num_atoms; $i++){
		$x[$i]=&to_fraction($x[$i]);
		$y[$i]=&to_fraction($y[$i]);
		$z[$i]=&to_fraction($z[$i]);
	}

#moto no kousi deno sisuu
#cartesian ni suru
	my (@r, @angle);
	for (my $i=0; $i<$num_atoms; $i++){
		($r[$i],$angle[$i])=&polar($xc[$i],$yc[$i],$zc[$i],@latvec_rotated[0..5]);
	}

#sort: angle sort
	my @order_key=&order_key_number(@angle);
	@x=@x[@order_key];
	@y=@y[@order_key];
	@z=@z[@order_key];
	@r=@r[@order_key];
	@angle=@angle[@order_key];
#chikai angle wo precise
	for (my $i=1; $i<=$#r; $i++){
		$angle[$i]=$angle[$i-1] if (&precise($angle[$i]-$angle[$i-1],0)==0);
	}

#sort: r sort
	@order_key=&order_key_number(@r);
	@x=@x[@order_key];
	@y=@y[@order_key];
	@z=@z[@order_key];
	@r=@r[@order_key];
	@angle=@angle[@order_key];
#saisyou angle wo 0, r wo precise
	my $baseangle=$angle[0];
	for (my $i=0; $i<$num_atoms; $i++){
		$angle[$i]-=$baseangle;
		$angle[$i]+=360 if ($angle[$i] < 0);
		$r[$i]=$r[$i-1] if (&precise($r[$i]-$r[$i-1],0)==0);
	}
#angle sort
	@order_key=&order_key_number(@angle);
	@x=@x[@order_key];
	@y=@y[@order_key];
	@z=@z[@order_key];
	@r=@r[@order_key];
	@angle=@angle[@order_key];

	for (my $i=$#r; $i>0; $i--){
#onaji angle nara ookii r wo korosu
		if (&precise($angle[$i]-$angle[$i-1],0,1E-6)==0){
			splice(@x,$i,1);
			splice(@y,$i,1);
			splice(@z,$i,1);
			splice(@r,$i,1);
			splice(@angle,$i,1);
		}
	}
#final r sort
	@order_key=&order_key_number(@r);
	@x=@x[@order_key];
	@y=@y[@order_key];
	@z=@z[@order_key];
	@r=@r[@order_key];
	@angle=@angle[@order_key];
#nagasa 12 made
	for (my $i=0; ($r[$i]<12 && $i < $#r); $i++){
		print ("$x[$i] $y[$i] $z[$i] : ");
		printf("%2.3f %3.1f \n", $r[$i], $angle[$i]);
	}
	exit;
}

sub get_terrace_vicinal_vector{
=pod
hkl slab, edge houkou hedge, kedge, ledge (he,ke,le) ga ataerareta baai, 
terace houkou hterrace, kterrace, lterrace (ht,kt,lt) to vicinal houkou
hvicinal, kvicinal, kvicinal (hv,kv,lv) no kouho wo kaesu
step ni tsuite ha 

edge houkou, terrace houkou ni tsuite ha
(h k l) * (he ke le) = 0
(h k l) * (ht kt lt) = 0
(h k l) * (hv kv lv) != 0
ga kyousei sareru

@_=(h,k,l,he,ke,le)
=cut

#hikisuu check
	my ($h, $k, $l)=@_[0..2];
	die ("H K L no H ga seisuu de nai: $h \n") if ($h !~ /^[-]?[0-9]+$/);
	die ("H K L no K ga seisuu de nai: $k \n") if ($k !~ /^[-]?[0-9]+$/);
	die ("H K L no L ga seisuu de nai: $l \n") if ($l !~ /^[-]?[0-9]+$/);
	die ("H K L ga zenbu 0 \n") if (&norm(&diff_vv($h,$k,$l,0,0,0)) == 0);
	die ("get_terrace_vicinal_vector: hikisuu ga 6 ka 7 yori chiisai: @_") if ($_[5] eq "");
	die ("get_terrace_vicinal_vector: hikisuu ga 6 ka 7 yori ookii: @_") if ($_[7] ne "");
	die ("get_terrace_vicinal_vector: hikisuu 7 ga asis de nai: @_") if (($_[6] ne "") && ($_[6] ne "asis"));

#conventional cell ni suru (hituyou nara)
	&get_bposcar if ($_[6] ne "asis");

#	my @hkle=&frac2dec_gcm(@_[3..5]);
	my @hkle=&frac2dec_array(@_[3..5]);
	my @hkle_integer=&frac2dec_gcm(@_[3..5]);
	die ("get_terrace_vicinal_vector: hkl(out-of-plane) @_[0..2] to hkl(edge) @_[3..5] no naiseki ga 0 de nai\n") if (&precise(&product_vv(@_[0..2],@hkle),0) != 0);
#hkle cartesian
	my @hkle_cart=&product_vm(@hkle,@latvec);


#zantei hv,kv,lv
	my @zhklv=&cross_vv(@hkle,@_[0..2]);
	my @latvec_first=@latvec;

#supercell;
#edge houkou to hkl houkou ha 2bai
#terrace houkou to hkl houkou ha $search_range
	my @zhklv_proj=&projection_normal_vector3D(@hkle_cart,&product_vm(@zhklv,@latvec));
	my $search_range=&ceiling(20/&norm(@zhklv_proj));

	&make_empty_cell(&product_vs(@zhklv,$search_range),&product_vs(@hkle,$h,$k,$l,2));

	my ($miny, $minz)=(999,999);
#empty cell no gensi suu
	my $num_atoms_bak=&round($num_atoms/4/$search_range);
	for (my $i=$num_atoms-1; $i>=0; $i--){
#miny on y axis sagashi
		$miny=$y[$i] if (($x[$i] == 0) && ($z[$i] == 0) && ($miny > $y[$i]) && ($y[$i] > 0));
#y zurashi
		$y[$i]-=0.5;
#minz sagashi
		$z[$i]=&precise($z[$i],$minz);
		$minz=$z[$i] if (($z[$i] > 0) && ($minz > $z[$i]));
	}

#iranai yz no gensi wo kesu
#moshikuha terrace site (z=0) to step site (z>0) no hutatsu ni wakery
	@namespecies=("T","S");
	@numspecies=(0,0);
	my (@x1, @x2, @y1, @y2, @z1, @z2);
	for (my $i=0; $i<$num_atoms; $i++){
		$y[$i]=&precise($y[$i],0.5*$miny);
		$y[$i]=&precise($y[$i],-0.5*$miny);
		$x[$i]=&precise($x[$i],0);
		$y[$i]=&precise($y[$i],0);
		$z[$i]=&precise($z[$i],0);

		if (($y[$i]>0.5*$miny) || ($y[$i] <=-0.5*$miny) || ($z[$i] > $minz)) {
# sakujyo suru gensi
		} elsif ($z[$i]==0){
			if (($x[$i] !=0) || ($y[$i] !=0)){
#terrace
				my @vec_cart=&product_vm($x[$i],$y[$i],$z[$i],@latvec);
				my $inner=&round(&product_vv(@vec_cart,@hkle_cart)/&product_vv(@hkle_cart,@hkle_cart));
				$y[$i]-=$inner*0.5;
				push (@x1,$x[$i]);
				push (@y1,$y[$i]);
				push (@z1,$z[$i]);
				$numspecies[0]++;
			}
		} else {
#step

			my @vec_cart=&product_vm($x[$i],$y[$i],$z[$i],@latvec);
			my $inner=&round(&product_vv(@vec_cart,@hkle_cart)/&product_vv(@hkle_cart,@hkle_cart));
			$y[$i]-=$inner*0.5;
			push (@x2,$x[$i]);
			push (@y2,$y[$i]);
			push (@z2,$z[$i]);
			$numspecies[1]++;
		}
	}
#merge
	@x=(@x1,@x2);
	@y=(@y1,@y2);
	@z=(@z1,@z2);
	$num_atoms=$numspecies[0]+$numspecies[1];

#plane normal
	my $height=&det3(@latvec)/&norm(&cross_vv(@latvec[0..5]))*$minz;
	printf ("step height: %2.2f\n",$height);
#lattice vector wo kaeru
	&d2c;
	@latvec=@latvec_first;
	my @xc=@x;
	my @yc=@y;
	my @zc=@z;
	my (@r_true, @r_proj, @angle, @len, @summary);
	&c2d;
	for (my $i=0; $i<$num_atoms; $i++){
		$summary[$i]=(&to_fraction($x[$i])." ".&to_fraction($y[$i])." ".&to_fraction($z[$i]));
		$r_true[$i]=&norm($xc[$i],$yc[$i],$zc[$i]);
		$r_proj[$i]=&norm(&projection_normal_vector3D(@hkle_cart,$xc[$i],$yc[$i],$zc[$i]));
	}
#sort,print
#terrace vector
	print ("terrace vector: \n");
	my @summary1=@summary[0..($numspecies[0]-1)];
	my @r_proj1=@r_proj[0..($numspecies[0]-1)];
	my @r_true1=@r_true[0..($numspecies[0]-1)];
	my @order_key=&order_key_number(@r_proj1);
	@summary1=@summary1[@order_key];
	@r_proj1=@r_proj1[@order_key];
	@r_true1=@r_true1[@order_key];
	for (my $i=0; $i<=$#r_proj1; $i++){
		print ("$summary1[$i] shaei_nagasa ");
		printf ("%3.2f true_nagasa %3.2f\n",$r_proj1[$i],$r_true1[$i]);
	}

#vicinal vector
	print ("vicinal vector: \n");
	my @summary1=@summary[$numspecies[0]..$num_atoms-1];
	@r_proj1=@r_proj[$numspecies[0]..$num_atoms-1];
	@r_true1=@r_true[$numspecies[0]..$num_atoms-1];
#angle
	my @proj_vector1=&projection_normal_vector3D(@hkle_cart,$xc[0],$yc[0],$zc[0]);
	for (my $i=$numspecies[0]; $i<$num_atoms;$i++){
		my @proj_vector2=&projection_normal_vector3D(@hkle_cart,$xc[$i],$yc[$i],$zc[$i]);
		$angle[$i-$numspecies[0]]=&angle_vv(0,0,0,@proj_vector1,@proj_vector2);
		$len[$i-$numspecies[0]]=&norm($xc[$i],$yc[$i],$zc[$i]);
	}


	my @x2=&product_vs(@x[$numspecies[0]..$num_atoms-1],$num_atoms_bak);
	my @y2=&product_vs(@y[$numspecies[0]..$num_atoms-1],$num_atoms_bak);
	my @z2=&product_vs(@z[$numspecies[0]..$num_atoms-1],$num_atoms_bak);
	@order_key=&order_key_number(@r_proj1);
	@summary1=@summary1[@order_key];
	@r_proj1=@r_proj1[@order_key];
	@r_true1=@r_true1[@order_key];
	@angle=@angle[@order_key];
	@x2=@x2[@order_key];
	@y2=@y2[@order_key];
	@z2=@z2[@order_key];
	@len=@len[@order_key];
	for (my $i=0; $i<=$#r_proj1; $i++){
		if ($angle[$i] < 90){
			print ("$summary1[$i] shaei_nagasa ");
			printf ("%3.2f true_nagasa %3.3f ",$r_proj1[$i],$r_true1[$i]);
			my @c=&cross_vv($x2[$i],$y2[$i],$z2[$i],@hkle_integer);
			@c=&round_array(&product_vs(@c,1/&gcmmany(@c)));
			@c=&minimum_array(@c);
			print ("surface @c angle(sigma) ");
			printf ("%3.1f sigma+tau=90_shaei_terrace_nagasa %3.2f\n",$angle[$i], $r_proj1[$i]*&cosd($angle[$i]));
		}
	}
	exit;
}



sub make_edge_slab{

=pod
edge ga aru slab no model wo tsukuru
edge vector: @_[0..2] terrace vector: @_[3..5] vicinal vector: @_[6..8]
slab thickness: $_[9] vacuum thickness: $_[10] asis: $_[11] no facet: $_[12] ("nofacet")
=cut
#kyousei d
	if ($dc eq "C"){
		&c2d;
		$dc="D";
	}
	&get_bposcar if ($_[11] ne "asis");
#get c-axis
	my @edge_dec=&frac2dec_array(@_[0..2]);
	my @terrace_dec=&frac2dec_array(@_[3..5]);
	my @vicinal_dec=&frac2dec_array(@_[6..8]);
	print ("decimal direct: edge @edge_dec terrace @terrace_dec vicinal @vicinal_dec \n");
	my @c=&minimum_array("keep_sign",&cross_vv(&frac2dec_gcm(@_[6..8]),&frac2dec_gcm(@_[0..2])));

#kezuru kakudo kakunin
	my @edge_vector_cartesian=&product_vm(@edge_dec,@latvec);
	my @vicinal_vector_cartesian=&product_vm(@vicinal_dec,@latvec);
	my @terrace_vector_cartesian=&product_vm(@terrace_dec,@latvec);
	my @step_vector_cartesian=&product_vm(&diff_vv(@vicinal_dec,@terrace_dec),@latvec);
#supercell
	&get_supercell(@vicinal_dec,@edge_dec,@c);

	my ($sigma, $tau);
#vector kakunin
	if ($_[12] ne "nofacet"){
		print ("cartesian: edge @edge_vector_cartesian vicinal @vicinal_vector_cartesian terrace @terrace_vector_cartesian step @step_vector_cartesian \n");
		my @vicinal_vector_cartesian_projected=&projection_normal_vector3D(@edge_vector_cartesian,@vicinal_vector_cartesian);
		my @terrace_vector_cartesian_projected=&projection_normal_vector3D(@edge_vector_cartesian,@terrace_vector_cartesian);
		my @step_vector_cartesian_projected=&projection_normal_vector3D(@edge_vector_cartesian,@step_vector_cartesian);
		my $norm_vicinal_vector_projected=&norm(@vicinal_vector_cartesian_projected);
		my $norm_terrace_vector_projected=&norm(@terrace_vector_cartesian_projected);
		my $norm_step_vector_projected=&norm(@step_vector_cartesian_projected);

		print ("projected_vicinal @vicinal_vector_cartesian_projected projected_terrace @terrace_vector_cartesian_projected projected_step @step_vector_cartesian_projected \n");
		printf ("lengths of projected: vicinal %2.2f terrace %2.2f step %2.2f \n",$norm_vicinal_vector_projected, $norm_terrace_vector_projected, $norm_step_vector_projected);
# vicinal^2 < terrace^2 + step^2 nara kanarazu 180-sigma-tau <= 90 deg
# vicinal^2 < terrace^2 + step^2 nara donkaku sankakkei de 0<sigma<90, 0<tau<90
		my $tolerance=1E-4;
		die ("make_edge_slab: projection: vicinal^2 < terrace^2 + step^2 \n") if ($norm_vicinal_vector_projected*$norm_vicinal_vector_projected < ($norm_terrace_vector_projected*$norm_terrace_vector_projected + $norm_step_vector_projected*$norm_step_vector_projected-$tolerance));

		$sigma=&angle_vv(0,0,0,@vicinal_vector_cartesian_projected,@terrace_vector_cartesian_projected);
		$tau=&angle_vv(0,0,0,@vicinal_vector_cartesian_projected,@step_vector_cartesian_projected);
		print ("sigma $sigma tau $tau \n");
		die ("make_edge_slab: not 0 < sigma < 90 : $sigma (denai hazu)\n") if (&within($sigma,0,90) == 0);
		die ("make_edge_slab: not 0 < tau < 90 : $tau (denai hazu) \n") if (&within($tau,0,90) == 0);
		die ("make_edge_slab: (sigma+tau) > 90 : sigma $sigma tau $tau (denai hazu) \n") if (($sigma+$tau)>90.001);
		print ("make supercell: @vicinal_dec @edge_dec @c \n");
	}

#c houkou ni chiisaku suru
	my @isometry=&get_isometry_translation;
	my @minz=(0,0,1);
	for (my $i=0; $i<$#isometry; $i+=3){
		@minz=@isometry[$i..($i+2)] if (($isometry[$i+2] > 0) && ($isometry[$i+2] < $minz[2]));
	}
	for (my $i=0; $i<3; $i++){
		$minz[$i]*=($#isometry+1);
		$minz[$i]=&round($minz[$i]);
		$minz[$i]/=($#isometry+1);
	}
#c jiku henkou
	&d2c;
	@latvec[6..8]=&product_vm(@minz,@latvec);
	&maxortho;
	&c2d;
	&unique();
#surface_isometry
	my ($numcenters,$shift)=&get_surface_slab_info;
	die ("Polar surface, koreijyou syori wo shimasen\n") if ($numcenters eq "Polar");
	my $slab_min_thick=$_[9];
	my $vac_min_thick=$_[10];
	my $vol=&det3(@latvec);
	my $prim_thick=$vol/&norm(&cross_vv(@latvec[0..5]));
	my $unit_thick=$prim_thick/$numcenters*2;
	my @latvec_backup=@latvec;
	my @x_backup=@x;
	my @y_backup=@y;
	my @z_backup=@z;
	my @w_backup=@w;
	my $num_atoms_backup=$num_atoms;
	my @numspecies_backup=@numspecies;

#slab no sozai wo tsukuru
#prim_thick and unit_thick in angstrom, shift in hkl-prim cell

	printf ("Saisyou slab thickness %2.2f Saisyou vacuum thickness %2.2f \n", $_[9],$_[10]);
	printf ("Primitive cell thickness %2.2f unit thickness %2.2f shift zPc %0.2f \n", $prim_thick, $unit_thick, $shift);

#half-integer=1, integer=2
	for (my $integer=1; $integer<=2; $integer++){
#center shift no=0, yes=1
		for (my $shiftcenter=0; $shiftcenter<=1; $shiftcenter++){
#slab model type (1-4)
			my $modeltype=($integer-1)*2+$shiftcenter;
#slab thickness (angstrom)
			my $slab_thick;
			if ($integer==1){
#half-integer
				$slab_thick=(&ceiling($slab_min_thick/$unit_thick-0.5)+0.5)*$unit_thick;
			} else {
#integer
				$slab_thick=&ceiling($slab_min_thick/$unit_thick)*$unit_thick;
			}
#hkl n-supercell
			my $n=&ceiling(($slab_thick+$vac_min_thick)/$prim_thick);
#cell thickness (angstrom)
			my $cell_thick=$n*$prim_thick;
			printf ("Slab thickness %2.2f(%0.2f) Cell thickness %2.2f N $n ",$slab_thick,$slab_thick/$cell_thick,$cell_thick);
#hkl-primitive cell de kangaeru
#center no ichi - tadasii atai (prim cell)
			my $realcenter=$shift+$shiftcenter/$numcenters;
#lower bound (prim cell)
			my $low_boundary=$realcenter-$integer/$numcenters/2;
			$low_boundary-=&floor($low_boundary);
			my $zminus=$low_boundary/$n;
			my $zplus=$zminus+$slab_thick/$cell_thick;
#makeslab
			printf ("z- %0.6f z+ %0.6f \n",$zminus,$zplus);

			&get_supercell(1,1,$n);

			&make_slab($zminus-1E-6,$zplus+1E-6);
#maxortho: use c
			&d2c;
			&maxortho;
			&c2d;
			&unique();

#saidai z ni chikai genshi (0 start)
#kokode osaete oku
			my @highz_atom;
			my $highz=&max(@z);
			for (my $i=0; $i<$num_atoms;$i++){
				push (@highz_atom, $i) if (&within($z[$i],$highz-(0.25/$n-1E-6), $highz+1E-6));
			}
			die ("make_edge_slab: highz ni chikai gensi ga nai : highz= $highz numatoms= $num_atoms (denai hazu)\n") if ($highz_atom[0] eq "");
			printf ("kiten no kazu ha %2d \n",$#highz_atom+1);

#nonpolar isometry check
			my @all_isometry=&get_isometry;
#a jiku houkou translation check
			my $counta;
			for (my $i=0;$i<$#all_isometry;$i+=12){
				my @t=@all_isometry[$i..($i+11)];
				$counta++ if (&norm(&diff_vv(@t[0..8],@t[10..11],1,0,0,0,1,0,0,0,1,0,0)) < 1E-5);
			}
# a jiku houkou translation w1 no saishou chi

			my @use_isometry=(0,0,0,0,0,0,0,0,0,999,999,999);
#type: 0 for -1, 1 for 2y,21y, 2 for m,b
			my $type=999;
			for (my $i=0;$i<$#all_isometry;$i+=12){
				my @t=@all_isometry[$i..($i+11)];
				if (&norm(&diff_vv(@t[0..8],-1,0,0,0,-1,0,0,0,-1)) < 1E-5){
#is -1
					if ($type > 0){
						@use_isometry=@t;
						$type=0;
					}
					@use_isometry=@t if ($t[9] <(1/$counta+0.0001));
				} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],-1,0,0,0,1,0,0,-1)) < 1E-5){
#is 2y, 21y
					if ($type > 1){
						@use_isometry=@t;
						$type=1;
					}
					@use_isometry=@t if (($type == 1) && ($t[9] <(1/$counta+0.0001)));

				} elsif  (&norm(&diff_vv(@t[0..1], @t[3..4], @t[6..8], 1,0,0,1,0,0,-1)) < 1E-5){
#is m, b
					if (&oddeven(2*(2*$t[9]+$t[2]*$t[11])*$counta)==0){;
						if ($type > 2){
							@use_isometry=@t;
							$type=2;
						}
						@use_isometry=@t if ((2*$t[9]+$t[2]*$t[11]) <(1/$counta+0.0001));
					}
				}
			}
print ("Isometry choice:@use_isometry use_isometry $counta counta $type type \n");
#slab model wo hozon
			my @latvec_slab_backup=@latvec;
			my @x_slab_backup=@x;
			my @y_slab_backup=@y;
			my @z_slab_backup=@z;
			my @w_slab_backup=@w;
			my $num_atoms_slab_backup=$num_atoms;
			my @numspecies_slab_backup=@numspecies;

#isometry ga aru, remove edge!
			if ($type != 999){
#onaji slab range no tooshi bangou 
				my $top_atom_count=0;
				if ($_[12] eq "nofacet"){
					my $filename="POSCAR.nofacet.".$modeltype;
					&writePOSCAR ($filename);
				} else {
#kiten goto ni shori
					for (my $i=0; $i<=$#highz_atom; $i++){
#edge nashi slab model ni modosu
						@latvec=@latvec_slab_backup;
						@x=@x_slab_backup;
						@y=@y_slab_backup;
						@z=@z_slab_backup;
						@w=@w_slab_backup;
						$num_atoms=$num_atoms_slab_backup;
						@numspecies=@numspecies_slab_backup;

#remove edge: kiten nokosu
						&remove_edge($x[$highz_atom[$i]],$z[$highz_atom[$i]],$sigma,$tau,@use_isometry);

						if ($#numspecies_backup != $#numspecies){
							printf ("reject: genso suu ga kawaru: edge nashi %d edge ari %d", $#numspecies_backup+1, $#numspecies+1) ;
						} elsif (&norm(&diff_vv(&minimum_array(@numspecies_backup), &minimum_array(@numspecies))) != 0){
#off-stoichiometry
							my @bak= &minimum_array(@numspecies_backup);
							my @new= &minimum_array(@numspecies);
							print ("off-stoichiometry: edge nashi @bak edge ari @new \n") ;
							my $file_edge="POSCAR.edge.".$modeltype.".".$top_atom_count.".offstoic";
							&writePOSCAR ($file_edge);
							$top_atom_count++;
						} else {
# stoichimetry OK
#nonpolar isometry check
							@all_isometry=&get_isometry;
#a jiku houkou translation check
							$counta=0;
							for (my $i=0;$i<$#all_isometry;$i+=12){
								my @t=@all_isometry[$i..($i+11)];
								$counta++ if (&norm(&diff_vv(@t,@use_isometry))<1E-4);
#isometry ga nakunaru = totemo mazui!!
							}
							die ("make_edge_slab: gensi kezuruto isometry ga ushinawareru: @use_isometry \n") if ($counta != 1);
																	my $file_edge="POSCAR.edge.".$modeltype.".".$top_atom_count.".vasp";
							&writePOSCAR ($file_edge);
							print ("stoichiometry: $file_edge \n") ;
							$top_atom_count++;
						}
#remove edge: kiten korosu
						@latvec=@latvec_slab_backup;
						@x=@x_slab_backup;
						@y=@y_slab_backup;
						@z=@z_slab_backup;
						@w=@w_slab_backup;
						$num_atoms=$num_atoms_slab_backup;
						@numspecies=@numspecies_slab_backup;

						&remove_edge($x[$highz_atom[$i]],$z[$highz_atom[$i]]-(1E-8),$sigma,$tau,@use_isometry);

						if ($#numspecies_backup != $#numspecies){
							printf ("reject: genso suu ga kawaru: edge nashi %d edge ari %d", $#numspecies_backup+1, $#numspecies+1) ;
						} elsif (&same_array(&minimum_array(@numspecies_backup), &minimum_array(@numspecies))){
#off-stoichiometry
							my @bak= &minimum_array(@numspecies_backup);
							my @new= &minimum_array(@numspecies);
							print ("off-stoichiometry: edge nashi @bak edge ari @new \n") ;
							my $file_edge="POSCAR.edge.".$modeltype.".".$top_atom_count.".offstoic";
							&writePOSCAR ($file_edge);
							$top_atom_count++;
						} else {
# stoichimetry OK
#nonpolar isometry check
							@all_isometry=&get_isometry;
#a jiku houkou translation check
							$counta=0;
							for (my $i=0;$i<$#all_isometry;$i+=12){
								my @t=@all_isometry[$i..($i+11)];
								$counta++ if (&norm(&diff_vv(@t,@use_isometry))<1E-4);
#isometry ga nakunaru = totemo mazui!!
							}
							die ("make_edge_slab: gensi kezuruto isometry ga ushinawareru: @use_isometry \n") if ($counta != 1);
																	my $file_edge="POSCAR.edge.".$modeltype.".".$top_atom_count.".vasp";
							&writePOSCAR ($file_edge);
							print ("stoichiometry: $file_edge \n") ;
							$top_atom_count++;
						}
					}
				}
			}
#termination wo kaeta node mdosu
			@latvec=@latvec_backup;
			@x=@x_backup;
			@y=@y_backup;
			@z=@z_backup;
			@w=@w_backup;
			$num_atoms=$num_atoms_backup;
			@numspecies=@numspecies_backup;
		}
	}

}

#edge wo kezuru
#hikisuu: kiten $x, $z, sigma, tau, @isometry (4-15)
sub remove_edge{
#epsilon, in angstrom
	my $eps=1E-4;
#position shift, angle
	my ($xshift, $zshift, $sigma, $tau)=@_[0..3];
#cartesian no dodai wo youi
	my @x_temp=@x;
	my @z_temp=@z;
	my @kill_atoms;
#cell geometry
	my $a=&norm(&projection_normal_vector3D(@latvec[3..5],@latvec[0..2]));
	my $c=&norm(&projection_normal_vector3D(@latvec[3..5],@latvec[6..8]));
	my $beta=&angle_vv(0,0,0,&projection_normal_vector3D(@latvec[3..5],@latvec[0..2]),&projection_normal_vector3D(@latvec[3..5],@latvec[6..8]));
	for (my $i=0; $i < $num_atoms; $i++){
		$x_temp[$i]-=$xshift;
		$z_temp[$i]-=$zshift;
		$x_temp[$i]++ if ($x_temp[$i] <0);
#cartesian
		my $x_cart=$a*$x_temp[$i]+$c*&cosd($beta)*$z_temp[$i];
		$x_cart=&precise($x_cart,0);
		while ($x_cart >= $a){
			$x_cart-=$a;
			$x_cart=&precise($x_cart,0);
		}
		while ($x_cart <0){
			$x_cart+=$a;
			$x_cart=&precise($x_cart,0);
		}
		my $z_cart=$c*&sind($beta)*$z_temp[$i];
#kezuru
#takasugi
		if ($z_cart > 0){
			push (@kill_atoms, $i+1);
			push (@kill_atoms, &isometry_image(@_[4..15],$i+1));
		} elsif ((&angle_vv(0,0,0,1,0,0,$x_cart-$eps,0,$z_cart)<$sigma) && (&angle_vv($a,0,0,0,0,0,$a-$x_cart-$eps,0,$z_cart)<$tau)){
#hutsuu ni kezuru
			push (@kill_atoms, $i+1);
			push (@kill_atoms, &isometry_image(@_[4..15],$i+1));
		}
	}
	print ("kill_atoms @kill_atoms \n");
	&remove_atoms(@kill_atoms);

}

#saikousei you no primitive cell wo tsukuru
sub get_primitive_for_surface_reconstruction{
#primitive
	my $tolerance=1E-8;
	&get_hkl_cell("primitive",@_);
	my @surface_isometry=&get_surface_isometry;
#saiyou suru isometry wo kakunin
#inversion (yuusen) mosikuha, mirror
	my @flag=(999,999,999,999);
	for (my $i=0;$i<$#surface_isometry;$i+=12){
		my @t=@surface_isometry[$i..($i+11)];
		my $test1=&oddeven(2*$t[9]+$t[2]*$t[11]);
		my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);

		if (&norm(&diff_vv(@t[0..8],-1,0,0,0,-1,0,0,0,-1)) < 1E-4){
		#-1
			print ("found -1! @t[9..11]\n");
			@t=&isometry_real(@t);
			@flag=(-1,@t[9..11]) if ($t[9] < $flag[1]);
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$test1,$test2,1,0,0,1,0,0,-1,0,0)) < 1E-4){
		#m
			print ("found m! $t[2] $t[5] @t[9..11]\n");
			@t=&isometry_real(@t);
			@flag=(0,@t[9..11]) if ($t[0] > 0);
		}
	}
	die ("inversion to mirror ga ryouhou nai \n") if ($flag[0] == 999);
	for (my $i=0; $i<$num_atoms; $i++){
		$x[$i]-=$flag[1]*0.5;
		$y[$i]-=$flag[2]*0.5;
		$z[$i]-=$flag[3]*0.5;
		$z[$i]++ if ($z[$i]<-$tolerance);
		$z[$i]=&precise($z[$i],0,$tolerance);
	}
	&d2c;
	&maxortho;
	&c2d;
	&site_in_cell;
	my @nonpolarity=&get_termination_polarity_nonpolarAB("");
	my $polartype=shift @nonpolarity;
	shift @nonpolarity;
	shift @nonpolarity;
	&maxortho;
	my $filename;
	if ($flag[0]==-1){
		$filename="POSCAR.INV.".$polartype.".";
	} else {
		$filename="POSCAR.MIR.".$polartype.".";
	}
	if ($polartype eq "A"){
		&writePOSCAR(${filename}."1.primitive");
		for (my $i=0; $i<$num_atoms; $i++){
			$z[$i]-=0.5;
			$z[$i]++ if ($z[$i]<-$tolerance);
		}
		&writePOSCAR(${filename}."2.primitive");
	} else {
		if (! grep { $_ == 0 } @nonpolarity){
			&writePOSCAR(${filename}."1.primitive");
		} else {
			print ("file 1 ga tsukurenai: nonpolar B, isometry no men ni genshi ga aru \n");
		}
		for (my $i=0; $i<$num_atoms; $i++){
			$z[$i]-=0.5;
			$z[$i]++ if ($z[$i]<-$tolerance);
		}
		if (! grep { $_ == 0.5 } @nonpolarity){
			&writePOSCAR(${filename}."2.primitive");
		} else {
			print ("file 2 ga tsukurenai: nonpolar B, isometry no men ni genshi ga aru \n");
		}
	}
}

sub get_termination_polarity{

=pod
termination wo kaeta toki no polarity wo kaesu
data wo temp_tsubo_polarity_data to shite kaesu
center no kazu, center no z houkou 0 kara no zure
return $flag,@termination,$thick,$numcenters,$shift, @return_list
saisyo ni $flag wo kaesu: X or (A,B,C)xn (example: AAAA, BC, CC)
$flag: A,B,C no kumiawase
@termination: 0,1,2,3,4 (nonpolar A,B: onaji termination ha onaji bangou
@return_list: $z1u, $z1c, low_boundary, polarity {Type_A, Type_B, Nonpolar Type_C , Non-stoichiometry}
=cut
	&get_hkl_cell("primitive",@_);
	my ($numcenters,$shift)=&get_surface_slab_info;
	die ("Polar surface, koreijyou syori wo shimasen\n") if ($numcenters eq "Polar");
	open OUTPUT, ">temp_tsubo_polarity_data";
	print OUTPUT ("(${_[0]}${_[1]}${_[2]})-primitive cell no jyouhou \n");
#kousi
	my @abcreal=&get_abc(@latvec);
	my @inplanelp=@abcreal[0..2];
	my $vol=&det3(@latvec);
	my $thick=$vol/&norm(&cross_vv(@latvec[0..5]));
	print OUTPUT ("Kousi no jyouhou\n");
	printf OUTPUT ("Real space a b gamma ((axb)*c)/|axb| : %2.2f %2.2f %2.2f %2.2f\n",@inplanelp[0..1],$abcreal[5],$thick); 
	my $det2d=&norm(&cross_vv(@latvec[0..5]))/6.283185307179586;
	printf OUTPUT ("Reciprocal space a* b*  : %2.2f %2.2f \n",&norm(@latvec[3..5])/$det2d,&norm(@latvec[0..2])/$det2d); 

#Model to termination no teigi
=pod
Model wo teigi suru ($i, $j)
I:   1,0
II:  1,1
III: 2,0
IV:  2,1
Gap (add $shift to all)
Tyouhuku termination no hanbetu ni siyou
G0: 0 
G1: 0 < z < 0.5/$numcenters
G2: 0.5/$numcenters
G3: 0.5/$numcenters < z < 1/$numcenters
G4: 1/$numcenters
Same if
I=III  : G3+G4 [A]
II=III : G2+G3 [B]
II=IV  : G0+G1 [C}
I=IV   : G1+G2 [D}

Nonpolar A, B de tyouhuku wo sakeru
Default:@termination=(1,2,3,4)
Nonpolar A, B de nai -> 0 ni suru
Onaji hyoumen nara onaji kigou ni suru

=cut
#gap no kakunin
	my @listz;
	my $tolerance=1E-8;
#	my $tolerance=1E-6;
#nonpolar termination: saidai yon shurui
#1,2,3,4 de zantei hyouki
#0: nonpolar A, B de wa nai
	my @termination=(1,2,3,4);
#G1 to G5, yes = 1, no = 0
	my @gap=(0,0,0,0,0);
	my $g0=$shift;
	my $g1=$shift+0.5/$numcenters;
	my $g2=$shift+1/$numcenters;
	for (my $m=0; $m<=$#z; $m++){
		$gap[0]=1 if (abs($z[$m]-$g0)<$tolerance);
		$gap[2]=1 if (abs($z[$m]-$g1)<$tolerance);
		$gap[4]=1 if (abs($z[$m]-$g2)<$tolerance);
		$gap[1]=1 if (&within($z[$m],$g0+$tolerance,$g1-$tolerance));
		$gap[3]=1 if (&within($z[$m],$g1+$tolerance,$g2-$tolerance));
	}
#termination check
#@atomexist: [A][B][C][D] ni gensi nasi = 0 ari = 1
	my @atomexist=(1,1,1,1);
	$atomexist[0]=0 if (&norm(@gap[3..4])==0);
	$atomexist[1]=0 if (&norm(@gap[2..3])==0);
	$atomexist[2]=0 if (&norm(@gap[0..1])==0);
	$atomexist[3]=0 if (&norm(@gap[1..2])==0);


#kokokara model wo kentou suru
#surface slab no boundary wo motomeru
	my @stoichiometry=&minimum_array(@numspecies);
#backup primtive
	my @latvec_orig=@latvec;
	my @x_orig=@x;
	my @y_orig=@y;
	my @z_orig=@z;
	my @w_orig=@w;
	my @numspecies_orig=@numspecies;
	my @namespecies_orig=@namespecies;

#1x1x3 supercell wo tsukuru
	&get_supercell(1,1,3);
#backup supercell
	my @latvec_backup=@latvec;
	my @x_backup=@x;
	my @y_backup=@y;
	my @z_backup=@z;
	my @w_backup=@w;
	my @numspecies_backup=@numspecies;
	my @return_list=($thick,$numcenters,$shift);
	my @output_list;
	print OUTPUT ("centers/unit (2/z_1u) $numcenters shift (z_Pc) $shift\n") ;
	print OUTPUT ("z_1u z_1c z_1- z_1+ Polarity TerminationID\n");

#$flag:A or B count
	my $flag;

#model wo tsukuru
	$tolerance/=3;
#haba
	for (my $i=1; $i<=2; $i++){
#takasa
		for (my $j=0; $j<2; $j++){
			my @polarity;
#slab no haba - tadasii atai 
			my $thickness=$i/$numcenters;
#center no ichi - tadasii atai 
			my $realcenter=$shift+$j/$numcenters;
#slab boundary	
			my @boundary=($realcenter-$thickness/2,$realcenter+$thickness/2);
#boundary wo ryougawa ireru (slab_long)
			my @slab_long=(($boundary[0]+1)/3-$tolerance,($boundary[1]+1)/3+$tolerance);
#boundary wo katahou ireru (slab_short)
			my @slab_short=(($boundary[0]+1)/3-$tolerance,($boundary[1]+1)/3-$tolerance);
#original slab thickness, slab center
			my $z1u=$boundary[1]-$boundary[0];
			my $z1c=($boundary[1]+$boundary[0])/2;
#ryoutan no gensi wo ireru (honrai no slab)
			&make_slab(@slab_long);
#zero gensi check
			my $zeroatom;
			for (my $m=0;$m<=$#numspecies_orig;$m++){
				$zeroatom++ if ($numspecies[$m] == 0);
			}
#zero gensi: slab wo atsuku suru (super long)
			if ($zeroatom > 0 ){
#slab tusukuri wo yarinaosu kara moto ni modosu
				@x=@x_backup;
				@y=@y_backup;
				@z=@z_backup;
				@w=@w_backup;
				@numspecies=@numspecies_backup;
				@namespecies=@namespecies_orig;

#aratamete  3-supercell kara slab wo kiridasu
				$thickness+=2/$numcenters;
				$boundary[1]+=2/$numcenters;
				my @slab_superlong=(($boundary[0]+1)/3-$tolerance,($boundary[1]+1)/3+$tolerance);
				&make_slab(@slab_superlong);
				@slab_short=(($boundary[0]+1)/3-$tolerance,($boundary[1]+1)/3-$tolerance);

			}
#stoichiometry = yes (=nonpolar A or B ga kakutei)
			if (&norm(&diff_vv(&minimum_array(@numspecies),@stoichiometry))==0){

#2679/11/1 subroutine ni suru
				my @test_polarity=&get_termination_polarity_nonpolarAB($flag);
				$flag=shift @test_polarity;
				@polarity=@test_polarity[0..1];
			} else {
#stoichiometry = no (=nonpolar C or non-stoichiometry ga kakutei)
#katahou no gensi wo hazusite stoichiometry check
				&make_slab(@slab_short);
				if (&norm(&diff_vv(&minimum_array(@numspecies),@stoichiometry))==0){
					$termination[($i-1)*2+$j]=($i-1)*2+$j+5;
					@polarity=("Nonpolar","Type_C");
					$flag.="C";
				} else {
					$termination[($i-1)*2+$j]=0;
					@polarity="Non-stoichiometry";
				}
			}
#non-polar nara data wo kaesu
			if ($#polarity != 0){
				my $lowboundary=$boundary[0]-&floor($boundary[0]);
				push (@return_list, $i, $j, $lowboundary, $polarity[1]);
			} else {
				$polarity[1]=" ";
			}
			@x=@x_backup;
			@y=@y_backup;
			@z=@z_backup;
			@w=@w_backup;
			@numspecies=@numspecies_backup;
			$num_atoms=&numatoms;
			push (@output_list, $z1u, $z1c, @boundary, @polarity);
		}
	}
#flag X check
	$flag="X" if ($flag eq "");
#termination kakunin 
#backup termination for nonpolar type C
	my @termination_backup=@termination;
=pod
Model same    $atomexist
I=III  :  [0]
II=III :  [1]
II=IV  :  [2]
I=IV   :  [3]
=cut

	if (($atomexist[0]+$atomexist[1]) == 0){
#no atom in [0],[1]
		$termination[2]=$termination[1]=$termination[0];
#I=II=III=IV if no atom in [2] or [3]
		$termination[3]=$termination[0] if (($atomexist[2]*$atomexist[3]) == 0);
	} elsif (($atomexist[2]+$atomexist[3]) == 0){
#no atom in [2],[3]
		$termination[3]=$termination[1]=$termination[0];
#I=II=III=IV if no atom in [0] or [1]
		$termination[2]=$termination[0] if (($atomexist[0]*$atomexist[1]) == 0);	} else {
#check IV then III
		$termination[3]=$termination[1] if ($atomexist[2] == 0);
		$termination[3]=$termination[0] if ($atomexist[3] == 0);
		$termination[2]=$termination[1] if ($atomexist[1] == 0);
		$termination[2]=$termination[0] if ($atomexist[0] == 0);
	}
#termination no bangou  wo tsumeru
=pod
1,2,3,4 de zantei hyouki
0: nonpolar A, B de wa nai
saisyo ha @termination=(1,2,3,4);
bangou ga 0 ni kawaru ka chiisaku natte iru
#sonzai sinai = 0 or 0 yori chiisai
=cut

	my @terminationexist=(0,0,0,0,0);
	for (my $i=0; $i<=3; $i++){
		$terminationexist[$termination[$i]]++;
	}
#termination 3 ga nai: then 4 -> 3
#termination 2 ga nai: then 4 -> 3, 3 -> 2
#termination 1 ga nai: then 4 -> 3, 3 -> 2, 2 -> 1
	for (my $i=3; $i>=1; $i--){
		if ($terminationexist[$i] == 0){
			for (my $j=0; $j<=3; $j++){
				$termination[$j]-- if ($termination[$j] > $i);
			}
		}
	}

#type C wo modosu
	for (my $i=0; $i<4; $i++){
		$termination[$i]=$termination_backup[$i] if ($termination_backup[$i] > 4);
	}

#print
	for (my $i=0; $i<4; $i++){
		printf OUTPUT ("%1.3f %1.3f %1.3f %1.3f ",@output_list[$i*6..($i*6+3)]);
		print OUTPUT (" $output_list[$i*6+4] $output_list[$i*6+5] ");
		print OUTPUT (" $termination[$i]") if ($termination[$i] != 0);
		print OUTPUT ("\n");
	}
	close OUTPUT;
#primitive ni modosu
	@x=@x_orig;
	@y=@y_orig;
	@z=@z_orig;
	@w=@w_orig;
	@latvec=@latvec_orig;
	@numspecies=@numspecies_orig;
	$num_atoms=&numatoms;

	return ($flag,@termination,@return_list);
}


#nonpolar A ka B wo kakunin
#$flag wo nyuuryoku
#A ka B ga tasareta $flag, "Nonpolar", "Type_A" ka "Type_B", @listz
sub get_termination_polarity_nonpolarAB{
	my @polarity="Nonpolar";
	my $flag=$_[0];
#Tasker 1 or 2 wo hantei
#z no list wo tsukuru
	my @listz=();
	for (my $m=0; $m<=$#z; $m++){
		my $newflag=0;
		for (my $n=0;$n<=$#listz;$n++){
			if ((&precise($z[$m],$listz[$n])) == $listz[$n]){
				$z[$m]=$listz[$n];
				$newflag++;
			}
		}
		$z[$m]=&precise($z[$m], 0, 1E-8);
		$z[$m]=&precise($z[$m], 0.5, 1E-8);
		push (@listz, $z[$m]) if ($newflag == 0);
	}
# tasker check 2d matrix: [z][genso]
	my @taskercheck;
	for (my $m=0; $m<=$#numspecies; $m++){
		for (my $k=0;$k<=$#listz;$k++){
			$taskercheck[$k][$m]=0;
		}
	}
	my $count=-1;
	for (my $m=0; $m<=$#numspecies; $m++){
		for (my $n=0; $n<$numspecies[$m]; $n++){
			$count++;
			for (my $k=0;$k<=$#listz;$k++){
				$taskercheck[$k][$m]++ if ($z[$count]==$listz[$k]);
			}
		}
	}
#hikaku
	$polarity[1]="Type_A";
	my @stoichiometry_base;
	for (my $m=0; $m<=$#numspecies; $m++){
		$stoichiometry_base[$m]=$taskercheck[0][$m];
	}
	@stoichiometry_base=&minimum_array(@stoichiometry_base);
	my @stoichiometry_check;
	for (my $k=1;$k<=$#listz;$k++){
		for (my $m=0; $m<=$#numspecies; $m++){
			$stoichiometry_check[$m]=$taskercheck[$k][$m];
		}
		@stoichiometry_check=&minimum_array(@stoichiometry_check);
		if (&same_array(@stoichiometry_check,@stoichiometry_base)){
			if ($polarity[1] eq "Type_A"){
				$polarity[1]="Type_B" ;
				$flag.="B";
			} 
		}
	}
	$flag.="A" if ($polarity[1] eq "Type_A");
	return ($flag, @polarity, @listz);
}


#Nonpolar type C slab no hyoumen no gensi wo hanbun kezuru
#cell wa hkl-primitive wo watasu
#hikisuu ha $n,$zminus,$zplus,$alternative (Y if alternative)
#layer group wo tsukau!
#assume D
#return 999 (tsukurenai)
#0 (tsukureta, alternative nashi)
#other 0~1 (tsukureta, alternative ari)
sub reconstruct_typeC_surface{
	my ($n, $zminus, $zplus, $alternative)=@_[0..3];
	my $bravais="other";

#nenno tame soseihi wo kakunin
	my @numspecies_bak=&minimum_array(@numspecies);

#2D Bravais lattice wo sagasu
	my @latpar=&get_abc(@latvec);
	if (abs($latpar[0]-$latpar[1])<1E-8){
#a=b
		if (&within($latpar[5],60-(1E-6),60+(1E-6))){
			&get_supercell(2,-1,0,0,1,0,0,0,1);
			$bravais="hp";
			print ("2D Bravais hp 2-101\n");
		} elsif (&within($latpar[5],120-(1E-6),120+(1E-6))){
			&get_supercell(2,1,0,0,1,0,0,0,1);
			$bravais="hp";
			print ("2D Bravais hp 2101\n");
		} elsif (&within($latpar[5],90-(1E-6),90+(1E-6))){
			$bravais="tp";
			print ("2D Bravais tp\n");
		} elsif (&within($latpar[5],60,90)){
			&get_supercell(1,1,0,-1,1,0,0,0,1);
			$bravais="oc";
			print ("2D Bravais oc 11-11\n");
		} elsif (&within($latpar[5],90,120)){
			&get_supercell(1,-1,0,1,1,0,0,0,1);
			$bravais="oc";
			print ("2D Bravais oc 1-111\n");
		} else {
			die ("reconstruct_typeC_surface: gamma_p ga okasii: $latpar[5] \n");
		}
	} else {
#a!=b
		my $t1=&product_vv(@latvec[0..5])*2;
		my $t2=&product_vv(@latvec[0..2],@latvec[0..2]);
		if (abs($t1+$t2)<1E-8){
			&get_supercell(-1,-2,0,1,0,0,0,0,1);
			$bravais="oc";
			print ("2D Bravais oc -1-201\n");
		} elsif (abs($t1-$t2)<1E-8){
			&get_supercell(1,-2,0,1,0,0,0,0,1);
			$bravais="oc";
			print ("2D Bravais oc 1-201\n");
		}
	}
#make slab
	&get_supercell(1,1,$n);	
	@latpar=&get_abc(@latvec);
	&make_slab($zminus-1E-6,$zplus+1E-6);
	my $wz=$zminus+$zplus;
#find center
	my $center;
	if ($bravais eq "hp"){
		print ("potential hp\n");
		$center=&nonpolarC_isometry_hp($wz);
	} elsif ($bravais eq "tp") {
		print ("potential tp\n");
		$center=&nonpolarC_isometry_tp($wz);
	} elsif ($bravais eq "oc") {
		print ("potential oc\n");
		$center=&nonpolarC_isometry_other($wz,"oc");
	} else {
		print ("other 2D Bravais lattice\n");
		$center=&nonpolarC_isometry_other($wz,"other");
	}
	return (999) if ($center == 999);
#top layer x_coordinate max gap
#max gap center: $x_gap[0]
	my @xlist;
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
		push (@xlist, $x[$i]) if (abs($z[$i]-$zplus) < 1E-6);
	}
	my @x_gap=&max_gap(@xlist);
	if ($alternative eq "Y"){
		if ($x_gap[2] ne ""){
			$x_gap[0]=$x_gap[2];
		} else {
			return (999);
		}
	}

#top center
	my $top_x_center=$x_gap[0]+0.25;
#0.25<=$top_x_center<0.75 no kukan ni gentei
	$top_x_center-=0.5 if ($top_x_center>=0.75);

#top kill range
	my @top_range=($top_x_center-0.25, $top_x_center+0.25);
	print ("top kill range : $top_x_center +/- 0.25 : $top_range[0] to $top_range[1] \n"); 
#bottom kill range 

#$bottom_x_center wo sagasu
	my $bottom_x_center;
	if ($center > 9){
		print ("reconstruct_typeC_surface: center $center > 9: nanika okasii !!!\n");
		return (999);
	} elsif ($center == 9){
		my @all_isometry=&get_isometry;
		for (my $i=0;$i<$#all_isometry;$i+=12){
			my @t=@all_isometry[$i..($i+11)];

			my $test1=&oddeven(2*(2*$t[9]+$t[2]*$t[11]));
			my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
			if (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$test2,1,0,0,1,0,0,-1,0)) < 1E-5){
#mirror, glide
				$center=$t[2] if ($t[9] < 0.5);
			} elsif (&norm(&diff_vv(@t[0..1],@t[3..8],$test1,1,0,0,-1,0,0,0,-1,0)) < 1E-5){
#2x
				$center=$t[2] if ($t[9] < 0.5);
			}
		}
# m, gb, 2a: xc'=xc+w13(zc-z-)
		print ("isometry type m/gb, w13 is $center  \n");
		$bottom_x_center=$top_x_center+$center*($zplus-$zminus)*0.5;		} else {
# -1 or 21x, xc'=-xc+w1
		print ("isometry type -1/2y/21y, w1 is $center  \n");
		$bottom_x_center=-1*$top_x_center+$center;
	}
	$bottom_x_center-=&floor($bottom_x_center);

#bottom layer kill range:
#0.25<$bottom_x_center<0.75: 0<kill range<1
#chigau baai(!)
#save range: ($bottom_x_center+0.25)-floor($bottom_x_center+0.25) to ($bottom_x_center-0.25)-floor($bottom_x_center-0.25)
	my $kill_flag=1;
	my @bottom_range;
	my $range1=$bottom_x_center-0.25;
	my $range2=$bottom_x_center+0.25;
	if (&within($bottom_x_center,0.25,0.75)){
		@bottom_range=($range1,$range2);
	} else {
		$kill_flag=0;
		@bottom_range=($range2-&floor($range2),$range1-&floor($range1));
	}
	print ("bottom range : @bottom_range \n"); 

#kill atom list
	my @remove_list;
	for (my $i=0; $i<$num_atoms; $i++){
#top layer kill
		if (abs($z[$i]-$zplus) < 1E-6){
			push (@remove_list, $i+1) if (&within($x[$i],@top_range));
		} elsif (abs($z[$i]-$zminus) < 1E-6){
			if ($kill_flag){
				push (@remove_list, $i+1) if (&within($x[$i],@bottom_range));
			} else {
				push (@remove_list, $i+1) if (! &within($x[$i],@bottom_range));
			}
		}
	}
#cell wo modosu

	&remove_atoms(@remove_list);


#nenno tame soseihi wo kakunin
	my @numspecies2=&minimum_array(@numspecies);
	if (&norm(&diff_vv(@numspecies_bak,@numspecies2)) != 0){
		die ("nonpolarC: soseihi ga chigau (arienai!) \n");
	}
#isometry wo kakunin
	my @all_isometry=&get_isometry;
	my $isometry_flag=0;
	print ("nonpolar wo tanpo suru isometry: ");
	for (my $i=0;$i<$#all_isometry;$i+=12){
		my @t=@all_isometry[$i..($i+11)];
		my $test1=&oddeven(2*$t[9]+$t[2]*$t[11]);
		my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
#inversion
		if (&norm(&diff_vv(@t[0..8],-1,0,0,0,-1,0,0,0,-1)) < 1E-5){
			print ("-1 [@t[9..10]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8], $test2,  -1,0,0,0,1,0,0,-1, 0)) < 1E-5){
#2y
			print ("2y [$t[9]] ");
			$isometry_flag=1;
		} elsif  (&norm(&diff_vv(@t[0..1],@t[3..8], $test1,  1,0,0,-1,0,0,0,-1, 0)) < 1E-5){
#2x
			print ("2x [$t[10]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8], $test2,  -1,0,0,0,1,0,0,-1, 1)) < 1E-5){
#21y
			print ("21y [$t[9]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..8], $test1,  1,0,0,-1,0,0,0,-1, 1)) < 1E-5){
#21x
			print ("21x [$t[9]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8], $test1,$test2,  1,0,0,1,0,0,-1, 0,0)) < 1E-5){
#m
			print ("m [@t[9..10]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8], $test1,$test2,  1,0,0,1,0,0,-1, 0,1)) < 1E-5){
#gb
			print ("gb [@t[9..10]] ");
			$isometry_flag=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8], $test1,$test2,  1,0,0,1,0,0,-1, 1,0)) < 1E-5){
#ga
			print ("ga [@t[9..10]] ");
			$isometry_flag=1;
		}
	}
	print ("\n");
	die ("reconstruct_typeC_surface: gensi kezutta ato top layer to bottom layer wo tsunagu isometry (primary isometry candidate) ga nai \n") if ($isometry_flag == 0);
	return ($x_gap[2]);

}

sub nonpolarC_isometry_hp{

#hex: yoko ni nagai double cell
	my @latvec_bak=@latvec;
	@latvec[0..2]=&product_vs(@latvec[0..2],0.58);
	my @all_isometry=&get_isometry;
#flags: -1, 2y, 2x, m
#data: wx@-1 [small wx], wy@-1 [small wy], wx @ 2y [small wx], wy @ 2x [small wy]
	my @flags=(0,0,0,0);
	my @data=(9,9,9,9);
	my $tolerance=1E-4;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		my @t=@all_isometry[$i..($i+11)];
		my $test1=&oddeven(2*$t[9]+$t[2]*$t[11]);
		my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
#inversion
		if (&norm(&diff_vv(@t[0..8],$t[11],-1,0,0,0,-1,0,0,0,-1,$t[11])) < 1E-5){
			print ("hp: found -1! vector @t[9..11]\n");
			$flags[0]=1;
			$data[0]=$t[9] if ($t[9]<$data[0]);
			$data[1]=$t[10] if ($t[10]<$data[1]);
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11], $test2,  -1,0,0,0,1,0,0,-1,$_[0], 0)) < $tolerance){
#2y
			print ("hp: found 2y! non-trivial matrix @t[2] vector  @t[9..11]\n");
			$flags[1]=1;
			$data[2]=$t[9] if ($t[9]<$data[2]);
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..8],$t[11], $test1,  1,0,0,-1,0,0,0,-1,$_[0], 0)) < $tolerance){
#2x
			print ("hp: found 2x! non-trivial matrix @t[5] vector  @t[9..11]\n");
			$flags[2]=1;
			$data[3]=$t[10] if ($t[10] < $data[3]);
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11], $test1,$test2,  1,0,0,1,0,0,-1,$_[0], 0,0)) < $tolerance){
#m
			print ("hp: found m! @t \n");
			$flags[3]=1;
		}
	}
	@latvec=@latvec_bak;
	print ("hp flags: @flags data: @data \n");
 	if (&same_array($flags[0],$flags[1],1,1) == 0){
#-1 + 2y
		print ("Principal Isometry 2y thru -1 \n");
		return ($data[2]);
	} elsif (&same_array($flags[0],$flags[2],1,1) == 0){
#-1 + 2x
		print ("Principal Isometry 2x thru -1 \n");
		&get_supercell(0,1,0,-1,0,0,0,0,1);
		return ($data[3]);
	} elsif ($flags[0]==1){
#-1
		print ("Principal Isometry -1 \n");
		return ($data[0]);
	} elsif ($flags[1]==1){
#2y
		print ("Principal Isometry 2y \n");
		return ($data[2]);
	} elsif ($flags[2]==1){
#2x
		print ("Principal Isometry 2x \n");
		&get_supercell(0,1,0,-1,0,0,0,0,1);
		return ($data[3]);
	} elsif ($flags[2]==1){
#m
		print ("Principal Isometry m \n");
		return (9);
	} else {
		print ("reconstruct_typeC_surface: potential hp 2D Bravais lattice, no principal isometry candidate \n") ;
		return (999);
	}
}


sub nonpolarC_isometry_tp{

#check flags: -1, 2y, 21y, 2+, 21+, 2-, 21-
#data: wx @ 2y, wx @ 21y
	my @all_isometry=&get_isometry;
	my @flags=(0,0,0,0,0,0,0);
	my @data=(9,9);
	my $tolerance=1e-4;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		my @t=@all_isometry[$i..($i+11)];
		my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
		my $test2p=&oddeven(2*($t[9]+$t[10]+$t[2]*$t[11]));
		my $test2m=&oddeven(2*($t[9]-$t[10]+$t[2]*$t[11]));
		print ("@t \n") if ($t[8] == -1);
#inversion
		if (&norm(&diff_vv(@t[0..8],$t[11],  -1,0,0,0,-1,0,0,0,-1,$_[0])) < $tolerance){
			print ("tp: found -1!  @t[9..10]\n");
			$flags[0]=1;
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11], $test2,  -1,0,0,0,1,0,0,-1,$_[0], 0)) < $tolerance){
#2y
			print ("tp: found 2y!  @t[9..10]\n");
			$flags[1]=1;
			$data[0]=$t[9];
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11], $test2,  -1,0,0,0,1,0,0,-1,$_[0], 1)) < $tolerance){
#21y
			print ("tp: found 21y!  @t[9..10]\n");
			$flags[2]=1;
			$data[1]=$t[9];
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11], $t[2],$test2p,  0,1,1,0,0,0,-1,$_[0], $t[5],0)) < $tolerance){
#2+
			print ("tp: found 2+\n");
			$flags[3]=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$t[2],$test2p,  0,1,1,0,0,0,-1,$_[0], $t[5],1)) < $tolerance){
#21+
			print ("tp: found 21+\n");
			$flags[4]=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$t[2],$test2m,  0,-1,-1,0,0,0,-1,$_[0], -$t[5],0)) < $tolerance){
#2-
			print ("tp: found 2-\n");
			$flags[5]=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$t[2],$test2m,  0,-1,-1,0,0,0,-1,$_[0], -$t[5],1)) < $tolerance){
#21-
			print ("tp: found 21-\n");
			$flags[6]=1;
		}

	}
	my $rotate=$flags[0]+$flags[3]+$flags[4]+$flags[5]+$flags[6];
#no rotation: no inverse, no diagonal
	if ($rotate == 0){
		if ($flags[1]==1){
			print ("2x1 supercell, principal isometry 2y \n");
			&get_supercell(2,1,1);
			return ($data[0]*0.5);
		} elsif ($flags[2]==1){
			print ("2x1 supercell, principal isometry 21y \n");
			&get_supercell(2,1,1);
			return ($data[1]*0.5);
		} else {
			print ("reconstruct_typeC_surface: potential tp 2D Bravais lattice ,no principal isometry candidate \n") ;
			return (999);
		}
	} else {
#45deg or 135deg rotate 
#kaiten kakudo wo kimeru
		if ($flags[3]==1){
#2+
			print ("45deg rotate \n");
			&get_supercell(1,1,0,-1,1,0,0,0,1);
		} elsif ($flags[4]==1) {
#2-
			print ("135deg rotate: hikakuteki mezurasii \n");
			&get_supercell(-1,1,0,-1,-1,0,0,0,1);
		} else {
#21+ or no 2-fold rotation
			print ("45deg rotate \n");
			&get_supercell(1,1,0,-1,1,0,0,0,1);
		}
#45deg or 135deg kaitenn ato
#flag: -1, 2y, 21y data: smallest wx @ -1, next wx @ -1, wx @ 2y, wx @ 21y
#2x, 21x ha kouryo sinakute ii
		@flags=(0,0,0);
		@data=(9,99,999,9999);
		@all_isometry=&get_isometry;
		for (my $i=0;$i<$#all_isometry;$i+=12){
			my @t=@all_isometry[$i..($i+11)];
			print ("@t \n") if ($t[8] == -1);
			my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
			if (&norm(&diff_vv(@t[0..8],$t[11],-1,0,0,0,-1,0,0,0,-1,$_[0])) < $tolerance){
#Inversion
				print ("tp_supercell: found -1!  @t[9..10]\n");
				$flags[0]=1;
				$data[1]=$t[9] if ($t[9] < $data[1]);
				@data[0..1]=($data[1], $data[0]) if ($data[1] < $data[0]);
			} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11],$test2,-1,0,0,0,1,0,0,-1,$_[0],0)) < $tolerance){
#2y
				print ("tp_supercell: found 2y!  $t[2] @t[9..10]\n");
				$flags[1]=1;
				$data[2]=$t[9] if ($t[9]<$data[2]);
			} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11],$test2,-1,0,0,0,1,0,0,-1,$_[0],1)) < $tolerance){
#21y
				print ("tp_supercell: found 21y! $t[2] @t[9..10]\n");
				$flags[2]=1;
				$data[3]=$t[9] if ($t[9]<$data[3]);
			}
		}
		print ("flags: @flags \n");
		print ("data: @data \n");
		if ($flags[0]==1){
			print ("Principal isometry -1 ");
			if ($data[2] == $data[0]){
				print ("penetrate 2y \n");
				return ($data[0]);
			} elsif ($data[2] == $data[1]){
				print ("penetrate 2y \n");
				return ($data[1]);
			} elsif ($data[3] == $data[0]){
				print ("penetrate 21y \n");
				return ($data[0]);
			} elsif ($data[3] == $data[1]){
				print ("penetrate 21y \n");
				return ($data[1]);
			} else {
				return ($data[0]);
			}
		} elsif ($flags[1]==1){
			print ("Principal isometry 2y\n");
			return ($data[2]);
		} elsif ($flags[2]==1){
			print ("Principal isometry 21y\n");
			return ($data[3]);
		} else {
			die ("reconstruct_typeC_surface: potential tp layer group, 45 or 135 deg rotate, no principal isometry candidate: totemo mazui! \n") ;
		}
	}
}


sub nonpolarC_isometry_other{
#flags: -1, 2y, 2x, 21y, 21x, m, ga, gb
#data= -1 wx(small), wy(small), wx(large), wy(large), 2y wx(small), wx(large), 2x wy(small), wy(large),21y wx(small), (large) ,21x wy(small), (large)
#already double (oc): $_[1] = "oc"
	my @all_isometry=&get_isometry;
	my @flags=(0,0,0,0,0,0,0,0);
	my @data=(9,9,99,99,8,88,8,88,8,88,8,88);
	my $tolerance=1E-4;
	for (my $i=0;$i<$#all_isometry;$i+=12){
		my @t=@all_isometry[$i..($i+11)];
		my $test1=&oddeven(2*$t[9]+$t[2]*$t[11]);
		my $test2=&oddeven(2*$t[10]+$t[5]*$t[11]);
		print ("@t \n") if ($t[8] == -1);
		#inversion 
		if (&norm(&diff_vv(@t[0..8],$t[11],-1,0,0,0,-1,0,0,0,-1,$_[0])) < $tolerance){
			print ("oc/op/mp: found -1! @t[9..10]\n");
			$flags[0]=1;
			if ($t[9] < $data[2]){
				if ($t[9] < $data[0]){
					$data[0]=$t[9];
				} else {
					$data[2]=$t[9];
				}
			}
			if ($t[10] < $data[3]){
				if ($t[10] < $data[1]){
					$data[1]=$t[10];
				} else {
					$data[3]=$t[10];
				}
			}
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11],$test2,-1,0,0,0,1,0,0,-1,$_[0],0)) < $tolerance){
		#2y
			print ("oc/op/mp: found 2y! $t[5] @t[9..10]\n");
			$flags[1]=1;
			if ($t[9] < $data[5]){
				if ($t[9] < $data[4]){
					$data[4]=$t[9];
				} else {
					$data[5]=$t[9];
				}
			}
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..8],$t[11],$test1,1,0,0,-1,0,0,0,-1,$_[0],0)) < $tolerance){
		#2x
			print ("oc/op/mp: found 2x! $t[2] @t[9..10]\n");
			$flags[2]=1;
			if ($t[10] < $data[7]){
				if ($t[10] < $data[6]){
					$data[6]=$t[10];
				} else {
					$data[7]=$t[10];
				}
			}
		} elsif (&norm(&diff_vv(@t[0..4],@t[6..8],$t[11],$test2,-1,0,0,0,1,0,0,-1,$_[0],1)) < $tolerance){
		#21y
			print ("oc/op/mp: found 21y! $t[5] @t[9..10]\n");
			$flags[3]=1;
			if ($t[9] < $data[9]){
				if ($t[9] < $data[8]){
					$data[8]=$t[9];
				} else {
					$data[9]=$t[9];
				}
			}
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..8],$t[11],$test1,1,0,0,-1,0,0,0,-1,$_[0],1)) < $tolerance){
		#21x
			print ("oc/op/mp: found 21x! $t[5] @t[9..10]\n");
			$flags[4]=1;
			if ($t[10] < $data[11]){
				if ($t[10] < $data[10]){
					$data[10]=$t[10];
				} else {
					$data[11]=$t[10];
				}
			}
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$test1,$test2,1,0,0,1,0,0,-1,$_[0],0,0)) < $tolerance){
		#m
			print ("op/mp: found m! $t[2] $t[5] @t[9..10]\n");
			$flags[5]=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$test1,$test2,1,0,0,1,0,0,-1,$_[0],1,0)) < $tolerance){
		#ga
			print ("op/mp: found ga! $t[2] $t[5] @t[9..10]\n");
			$flags[6]=1;
		} elsif (&norm(&diff_vv(@t[0..1],@t[3..4],@t[6..8],$t[11],$test1,$test2,1,0,0,1,0,0,-1,$_[0],0,1)) < $tolerance){
		#gb
			print ("op/mp: found gb! $t[2] $t[5] @t[9..10]\n");
			$flags[7]=1;
		}
	}
#flags: -1, 2y, 2x, 21y, 21x, m, ga, gb
#data= -1 wx(small), wy(small), wx(large), wy(large), 2y wx(small), wx(large), 2x wy(small), wy(large),21y wx(small), (large) ,21x wy(small), (large)
	print ("oc/op/mp flags @flags \n");
	print ("data @data \n");
#inversion
	if ($flags[0] == 1){
#thru 2y?
		if ($flags[1] == 1){
			if (($data[0] == $data[4])||($data[0] == $data[5])){
				print ("principal isometry 2y thru -1 [small wx]\n");
				if ($_[1] eq "other"){
				&get_supercell(2,1,1);
					return ($data[0]*0.5);
				} else {
					return ($data[0]);
				}
			} elsif (($data[2] == $data[4])||($data[2] == $data[5])){
				print ("principal isometry 2y thru -1 [large wx]\n");
				if ($_[1] eq "other"){
				&get_supercell(2,1,1);
					return ($data[2]*0.5);
				} else {
					return ($data[2]);
				}
			} 
		}
#thru 2x?
		if ($flags[2] == 1){
			if (($data[1] == $data[6])||($data[1] == $data[7])){
				print ("principal isometry 2x thru -1 [small wy] -> 2y thru -1\n");
				if ($_[1] eq "other"){
					&get_supercell(0,2,0,-1,0,0,0,0,1) ;
					return ($data[1]*0.5);
				} else	{
					&get_supercell(0,1,0,-1,0,0,0,0,1) ;
					return ($data[1]);
				}
			} elsif (($data[3] == $data[6])||($data[3] == $data[7])){
				print ("principal isometry 2x thru -1 [large wy] -> 2y thru -1\n");
				if ($_[1] eq "other"){
					&get_supercell(0,2,0,-1,0,0,0,0,1) ;
					return ($data[3]*0.5);
				} else	{
					&get_supercell(0,1,0,-1,0,0,0,0,1) ;
					return ($data[3]);
				}
			}
		}
#thru 2y1?
		if ($flags[3] == 1){
			if (($data[0] == $data[8])||($data[0] == $data[9])){
				print ("principal isometry 21y thru -1 [small wx]\n");
				if ($_[1] eq "other"){
					&get_supercell(2,1,1);
					return ($data[0]*0.5);
				} else {
					return ($data[0]);
				}
			} elsif (($data[2] == $data[8])||($data[2] == $data[9])){
				print ("principal isometry 21y thru -1 [large wx]\n");
				if ($_[1] eq "other"){
					&get_supercell(2,1,1);
					return ($data[2]*0.5);
				} else {
					return ($data[2]);
				}
			}
		}
#thru 2x1?
		if ($flags[4] == 1){
			if (($data[1] == $data[10])||($data[1] == $data[11])){
				print ("principal isometry 21x thru -1 [small wy] -> 21y thru -1\n");
				if ($_[1] eq "other"){
					&get_supercell(0,2,0,-1,0,0,0,0,1) ;
					return ($data[1]*0.5);
				} else	{
					&get_supercell(0,1,0,-1,0,0,0,0,1) ;
					return ($data[1]);
				}
			} elsif (($data[3] == $data[10])||($data[3] == $data[11])){
				print ("principal isometry 21x thru -1 [large wy] -> 21y thru -1\n");
				if ($_[1] eq "other"){
					&get_supercell(0,2,0,-1,0,0,0,0,1) ;
					return ($data[3]*0.5);
				} else	{
					&get_supercell(0,1,0,-1,0,0,0,0,1) ;
					return ($data[3]);
				}
			}
		}
		print ("principal isometry -1, no 2 nor 21 penetration \n");
		if ($_[1] eq "other"){
			&get_supercell(2,1,1); 
			return ($data[0]*0.5);
		} else {
			return ($data[0]);
		}
	} elsif ($flags[1] == 1){
		print ("principal isometry 2y\n");
		if ($_[1] eq "other"){
			&get_supercell(2,1,1);
			return ($data[4]*0.5);
		} else {
			return ($data[4]);
		}
	} elsif ($flags[2] == 1){
		print ("principal isometry 2x -> 2y\n");
		if ($_[1] eq "other"){
			&get_supercell(0,2,0,-1,0,0,0,0,1) ;
			return ($data[6]*0.5);
		} else	{
			&get_supercell(0,1,0,-1,0,0,0,0,1) ;
			return ($data[6]);
		}
	} elsif ($flags[3] == 1){
		print ("principal isometry 21y\n");
		if ($_[1] eq "other"){;
			&get_supercell(2,1,1);
			return ($data[8]*0.5);
		} else {
			return ($data[8]);
		}
	} elsif ($flags[4] == 1){
		print ("principal isometry 21x -> 21y\n");
		if ($_[1] eq "other"){
			&get_supercell(0,2,0,-1,0,0,0,0,1) ;
			return ($data[10]*0.5);
		} else	{
			&get_supercell(0,1,0,-1,0,0,0,0,1) ;
			return ($data[10]);
		}
	} elsif  ($flags[5] == 1){
#m
		&get_supercell(2,1,1) if ($_[1] eq "");
		print ("principal isometry m\n");
		return (9);
	} elsif ($flags[6] == 1){
#ga
		if ($_[1] eq "other"){
			&get_supercell(0,2,0,-1,0,0,0,0,1) ;
		} else	{
			&get_supercell(0,1,0,-1,0,0,0,0,1) ;
		}
		print ("principal isometry a -> b\n");
		return (9);
	} elsif ($flags[7] == 1){
#gb
		&get_supercell(2,1,1) if ($_[1] eq "other");
		print ("principal isometry b\n");
		return (9);
	} else {
		print ("reconstruct_typeC_surface: oc/mp/ap layer group, gensi kezutta ato top layer to bottom layer wo tsunagu isometry ga nai \n");
		return (999);
	}
}


sub facet_orientation{
#facet search
#arg: H K L maxH maxK maxL sgnumber(facet nomi hyouji no baai) library (energy wo dasu baai)
	my @tmatrix=&hkl_transformation_matrix(@_[0..2], "nonCMS");
	my @lattemp=&product_mm(@tmatrix[0..8],@latvec);
	my @a=&unit_v(@lattemp[0..2]);
	my @c=&unit_v(&cross_vv(@lattemp[0..5]));
	my @b=&unit_v(&cross_vv(@c,@a));
	my (@hh1, @kk1, @ll1, @theta1, @phi1);
	my (@hh2, @kk2, @ll2, @theta2, @phi2);
	my ($i, $j)=(-1, -1);
#h k l shirami tsubushi
	for (my $h=-$_[3]; $h<=$_[3]; $h++){
		for (my $k=-$_[4]; $k<=$_[4]; $k++){
			for (my $l=-$_[5]; $l<=$_[5]; $l++){
#check if no 0 0 0, nor not coprime
				my $flag=abs(&gcmmany($h,$k,$l));
#hkl ga @arg[0..2] to onaji nara 999
				$flag=999 if (&norm(&cross_vv(@_[0..2],$h,$k,$l))==0);
#flag = 1: OK, soreigai ha hajiku 
				if ($flag == 1){
					my @testmatrix=&hkl_transformation_matrix($h,$k,$l, "nonCMS");
					my @testlattemp=&product_mm(@testmatrix[0..8],@latvec);
					my @testc=&unit_v(&cross_vv(@testlattemp[0..5]));

					my $theta=&acos(&product_vv(@c,@testc));
#gyaku houkou no hyoumen wo hajiku
					if ($theta < 90){
						my @inplane;
						my @t=&polar(&projection_normal_vector3D(@c,@testc),@a,@b);
						my $phi=&round($t[1]*10);

#>=180 to <180 wo betu no retu de hyouji suru node hairetu wo kaeru
						if ($phi < 1800){
							$i++;
							$hh1[$i]=$h;
							$kk1[$i]=$k;
							$ll1[$i]=$l;
							$theta1[$i]=$theta;
							$phi1[$i]=$phi;
						} else {
							$j++;
							$hh2[$j]=$h;
							$kk2[$j]=$k;
							$ll2[$j]=$l;
							$theta2[$j]=$theta;
							$phi2[$j]=$phi-1800;
						}
					}
				}
			}
		}
	}
#zenbu hyouji mode
	if ($_[6] eq ""){
	#sort
		my @order_key=&order_key_number(@theta1);
		@hh1=@hh1[@order_key];
		@kk1=@kk1[@order_key];
		@ll1=@ll1[@order_key];
		@theta1=@theta1[@order_key];
		@phi1=@phi1[@order_key];
		@order_key=&order_key_number(@phi1);
		@hh1=@hh1[@order_key];
		@kk1=@kk1[@order_key];
		@ll1=@ll1[@order_key];
		@theta1=@theta1[@order_key];
		@phi1=@phi1[@order_key];
		@order_key=&order_key_number(@theta2);
		@hh2=@hh2[@order_key];
		@kk2=@kk2[@order_key];
		@ll2=@ll2[@order_key];
		@theta2=@theta2[@order_key];
		@phi2=@phi2[@order_key];
		@order_key=&order_key_number(@phi2);
		@hh2=@hh2[@order_key];
		@kk2=@kk2[@order_key];
		@ll2=@ll2[@order_key];
		@theta2=@theta2[@order_key];
		@phi2=@phi2[@order_key];

	#$i ga $j yori nagaku suru
		if ($j > $i){
			my @s=@hh2;
			@hh2=@hh1;
			@hh1=@s;
			@s=@kk2;
			@kk2=@kk1;
			@kk1=@s;
			@s=@ll2;
			@ll2=@ll1;
			@ll1=@s;
			@s=@theta2;
			@theta2=@theta1;
			@theta1=@s;
			@s=@phi2;
			@phi2=@phi1;
			@phi1=@s;
			($i,$j)=($j,$i);
		}
	#print
		print (" h k l theta phi h k l theta phi \n");
		for (my $x=0; $x<=$j; $x++){
			printf (" %2d %2d %2d %2.1f %3.1f  %2d %2d %2d %2.1f %3.1f (+180)\n",$hh1[$x], $kk1[$x], $ll1[$x], $theta1[$x], $phi1[$x]*0.1, $hh2[$x], $kk2[$x], $ll2[$x], $theta2[$x], $phi2[$x]*0.1);
		}
		for (my $x=$j+1; $x<=$i; $x++){
			printf (" %2d %2d %2d %2.1f %3.1f \n",$hh1[$x], $kk1[$x], $ll1[$x], $theta1[$x], $phi1[$x]*0.1);
		}
		exit;
	}
#facet only mode
#kyoutsuu suru phi
	my @unique_phi1=&unique_array(@phi1);
	my @unique_phi2=&unique_array(@phi2);
	my @orient1=@theta1;
	my @orient2=@theta2;
	my @orient_orig1=@theta1;
	my @orient_orig2=@theta2;
	my %temp_phi;
	my %energylibrary;
	my $surfE=99999;
	my $facet;
	my ($theta11, $theta22);
	$temp_phi{$_}++ for (@unique_phi1,@unique_phi2);
	my @share_phi = grep{ $temp_phi{$_} >= 2 } keys %temp_phi;
	for (my $x=$#phi1; $x>=0; $x--){
		my $save_flag=1;
		if (! grep { $_ eq $phi1[$x] } @share_phi){
			$save_flag=0;
		} else {
			my @temp=&to_unique_nonpolar($_[6],$hh1[$x],$kk1[$x],$ll1[$x]);
			if ($temp[0] eq "Nonpolar"){
				$save_flag=0;
			} else {
				$orient1[$x]=$temp[0]."".$temp[1]."".$temp[2];
				$orient_orig1[$x]=$hh1[$x]."".$kk1[$x]."".$ll1[$x];
			}
		}
		if ($save_flag == 0){
			splice(@hh1, $x, 1);
			splice(@kk1, $x, 1);
			splice(@ll1, $x, 1);
			splice(@theta1, $x, 1);
			splice(@phi1, $x, 1);
			splice(@orient1, $x, 1);
			splice(@orient_orig1, $x, 1);
		}
	}
	for (my $x=$#phi2; $x>=0; $x--){
		my $save_flag=1;
		if (! grep { $_ eq $phi2[$x] } @share_phi){
			$save_flag=0;
		} else {
			my @temp=&to_unique_nonpolar($_[6],$hh2[$x],$kk2[$x],$ll2[$x]);
			if ($temp[0] eq "Nonpolar"){
				$save_flag=0;
			} else {
				$orient2[$x]=$temp[0]."".$temp[1]."".$temp[2];
				$orient_orig2[$x]=$hh2[$x]."".$kk2[$x]."".$ll2[$x];
			}
		}
		if ($save_flag == 0){
			splice(@hh2, $x, 1);
			splice(@kk2, $x, 1);
			splice(@ll2, $x, 1);
			splice(@theta2, $x, 1);
			splice(@phi2, $x, 1);
			splice(@orient2, $x, 1);
			splice(@orient_orig2, $x, 1);
		}
	}
#merge
	for (my $x=$#hh2; $x>=0; $x--){
		push (@hh1, shift(@hh2));
		push (@kk1, shift(@kk2));
		push (@ll1, shift(@ll2));
		push (@theta1, 180-shift(@theta2));
		push (@phi1, shift(@phi2));
		push (@orient1, shift(@orient2));
		push (@orient_orig1, shift(@orient_orig2));
	}
#sort
	my @order_key=&order_key_number(@theta1);
	@hh1=@hh1[@order_key];
	@kk1=@kk1[@order_key];
	@ll1=@ll1[@order_key];
	@theta1=@theta1[@order_key];
	@phi1=@phi1[@order_key];
	@orient1=@orient1[@order_key];
	@orient_orig1=@orient_orig1[@order_key];
	@order_key=&order_key_number(@phi1);
	@hh1=@hh1[@order_key];
	@kk1=@kk1[@order_key];
	@ll1=@ll1[@order_key];
	@theta1=@theta1[@order_key];
	@phi1=@phi1[@order_key];
	@orient1=@orient1[@order_key];
	@orient_orig1=@orient_orig1[@order_key];
#energy mode de nai baai print
	if ($_[7] eq ""){
		print ("h k l unique theta phi \n");
		for (my $x=0; $x<=$#hh1; $x++){
			printf (" %2d %2d %2d $orient1[$x] %2.12f %3.1f\n",$hh1[$x], $kk1[$x], $ll1[$x], $theta1[$x], $phi1[$x]*0.1);
		}
		print ("Step kumiawase\n");
		print (" h1 k1 l1 unique1 theta1 h2 k2 l2 unique2 theta2 theta1+theta2 phi1\n");
	} else {
#library wo tsukuru
		open A, $_[7];
		while (my $t=<A>){
			my @list=&splitall($t);
			if (($list[0] ne "") && ($list[1] ne "")){
				$energylibrary{$list[0]}=$list[1];
			}
		}
		close A;
	}
#share_phi list
	@order_key=&order_key_number(@share_phi);
	@share_phi=@share_phi[@order_key];
	my @same_phi;
#for given phi
	for (my $x=0; $x<=$#share_phi; $x++){
		@same_phi = grep { $phi1[$_] eq $share_phi[$x] } 0 .. $#phi1
;
		my (@phi11, @phi22);
#theta <90 to >90 ni wakeru
		for (my $y=0; $y<=$#same_phi; $y++){
			if ($theta1[$same_phi[$y]] < 90){
				push (@phi11, $same_phi[$y]);
			} else {
				push (@phi22, $same_phi[$y]);
			}
		}
#kumiawase
		for (my $y=0; $y<=$#phi11; $y++){
			for (my $z=$#phi22; $z>=0; $z--){
				if (($theta1[$phi22[$z]]-$theta1[$phi11[$y]]) >= 90){
#print mode
					if ($_[7] eq ""){
						printf (" %2d %2d %2d $orient1[$phi11[$y]] %2.12f  ",$hh1[$phi11[$y]], $kk1[$phi11[$y]], $ll1[$phi11[$y]], $theta1[$phi11[$y]], );
						printf ("%2d %2d %2d $orient1[$phi22[$z]] %2.12f  %2.12f %3.1f\n", $hh1[$phi22[$z]], $kk1[$phi22[$z]], $ll1[$phi22[$z]], 180-$theta1[$phi22[$z]], 180-$theta1[$phi22[$z]]+$theta1[$phi11[$y]], $share_phi[$x]*0.1);
					} else {
#facet surfE
						if (exists($energylibrary{$orient1[$phi11[$y]]}) && exists($energylibrary{$orient1[$phi22[$z]]})){

							my $facetE=($energylibrary{$orient1[$phi11[$y]]}*&sind($theta1[$phi22[$z]])+$energylibrary{$orient1[$phi22[$z]]}*&sind($theta1[$phi11[$y]]))/&sind($theta1[$phi22[$z]]-$theta1[$phi11[$y]]);
							if ($facetE < $surfE){
								$surfE=$facetE;
								$facet=$orient1[$phi11[$y]]."/".$orient1[$phi22[$z]];
								$facet_orig=$orient_orig1[$phi11[$y]]."/".$orient_orig1[$phi22[$z]];
								$theta11=$theta1[$phi11[$y]];
								$theta22=$theta1[$phi22[$z]];
							}
						}
					}
				}
			}
		}
	}
#facet E hikaku mode only
	if ($_[7] ne ""){
		die ("Original @_[0..2] : facet no energy ga dasenai\n") if ($surfE == 99999);
		print ("Minimum_facet original $facet_orig unique $facet $surfE No_facet Original $_[0]$_[1]$_[2] Unique ");
		my ($h, $k, $l)=&to_unique_nonpolar($_[6],@_[0..2]);
		my $orient=${h}.${k}.${l};
		print ("$orient ");
		if (exists($energylibrary{$orient})){
			my $origE=$energylibrary{$orient};
			my $delta=$surfE-$origE;
			print ("$origE Delta $delta ");
			printf ("theta1 %3.2f theta2 %3.2f edgeangle %3.2f ",$theta11,180-$theta22, $theta22-$theta11);
		}
		print ("\n");
		exit;
	}
}




############
#
#  Band path
#
############


sub get_kpath_kpf_print{
#hikisuu: kou taishou ten toosi bangou (1 start), kou taishou ten bangou 
#(0 = Gamma), type (-, |, X no izureka), kou taisyou ten hairetu
# printf OUT ga tsukaeru koto!
	printf OUT (" % 2.7f % 2.7f % 2.7f K.%d %s %s\n", @_[$_[1]*4+3..$_[1]*4+5], $_[0], $_[$_[1]*4+6], $_[2]);
}

sub get_kpath_kpoints{
#hikisuu: $_[0] = kpoint no kankaku (units A-1), $_[1] = Kessyougaku or SC
#$_[2] = force time reversal symmetry (Y/N)
#kpf data wo tsubo_temp_get_kpath kara yomu
#KPOINT keisiki data wo aratamete tsubo_temp_get_kpath ni syuturyoku

#kantann na hikisuu check: 
	die ("Kpoint kankaku $_[0] ga kazu de nai\n") if ($_[0] =~ /[a-zA-Z]/);
	die ("Kpoint kankaku $_[0] ga 0 ika \n") if ($_[0] <= 0);
#kpf
	if ($_[1] eq "SC") {
		&get_kpath_kpf_SC;
	} else {
		&get_kpath_kpf($_[2]);
	}
	my @points;
	open KPT, "tsubo_temp_get_kpath";
	my $a=<KPT>;
	my $numkpoints=0;
	while ($a=<KPT>){
		chomp $a;
		my @b = split(/ +/,$a);
		$points[$numkpoints][0]=$b[6];
		$points[$numkpoints][1]=$b[1];
		$points[$numkpoints][2]=$b[2];
		$points[$numkpoints][3]=$b[3];
		$points[$numkpoints][8]=$b[4];
		$points[$numkpoints][9]=$b[5];
		$numkpoints++;
	}
	close KPT;
	
#@points matrix wo tsukuru
# Main matrix $points[site][item]
# item: 0: -|X 1-3: Fractional 4-6: Cartesian coord. 
# 7: Mae no ten kara no kyori 8: K.X 9: Symbol
	my @recilatvec=&reci_latvec(@latvec);
	my $sumlength=0;
	for (my $i=0;$i<$numkpoints;$i++){
		($points[$i][4],$points[$i][5],$points[$i][6])=
		    &product_vm($points[$i][1],$points[$i][2],$points[$i][3],@recilatvec);
		if ($i==0){
			$points[0][7]=0;
		} else {
			$points[$i][7]=sqrt(($points[$i][4]-$points[$i-1][4])**2+($points[$i][5]-$points[$i-1][5])**2+($points[$i][6]-$points[$i-1][6])**2);
			$sumlength+=$points[$i][7];
		}
	}

#syuturyoku 
	my @vicinal;
	my $writeflag=1;
	for (my $i=1;$i<$numkpoints;$i++){
#get number of points in segment
		$vicinal[$i]=int($points[$i][7]/$_[0]+0.5);
	}
	open OUT, "> tsubo_temp_get_kpath";
#saisyo no site
	printf OUT (" % 2.7f % 2.7f % 2.7f ", $points[0][1], $points[0][2], $points[0][3]);	
	print OUT ("0 ! $points[0][8] $points[0][9] $points[0][0] \n");
	for (my $i=1;$i<$numkpoints;$i++){
		if ($writeflag == 1){
			for (my $j=1;$j<$vicinal[$i];$j++){
				for (my $k=1; $k<=3;$k++){
					printf OUT (" % 2.7f",$points[$i][$k]*($j/$vicinal[$i])+$points[$i-1][$k]*(1-$j/$vicinal[$i]));
				}
				print OUT (" 0 \n");
			}
		}
		printf OUT (" %2.7f %2.7f %2.7f ", $points[$i][1], $points[$i][2], $points[$i][3]);	
		print OUT ("0 ! $points[$i][8] $points[$i][9] $points[$i][0] \n");
		$writeflag=0 if ($points[$i][0] eq "|");
		$writeflag=1 if ($points[$i][0] eq "-");
	}
	close OUT;
} 

sub get_kpath_phonopy{
#PHONOPY keisiki data wo aratamete tsubo_temp_get_kpath ni syuturyoku
#Kessyougaku to SC wo $_[0] de ataeru
#start
	my $band="BAND= ";
	my $band_labels="BAND_LABELS= ";
	my $gamma="\Gamma ";
	if ($_[0] eq "SC") {
		&get_kpath_kpf_SC;
	} else {
		&get_kpath_kpf;
	}
# kpf format kpoints wo yomu
	open KPT, "tsubo_temp_get_kpath";
	my $a=<KPT>;
	my $numkpoints=0;
	while ($a=<KPT>){
		chomp $a;
#$a format ha
#   0.0000000  0.0000000  0.0000000 K.1 Gamma -
		my @b = split(/ +/,$a);
		for (my $i=1;$i<=3;$i++){
			$b[$i]=&precise($b[$i],0);
		}
		$band.=$b[1]." ".$b[2]." ".$b[3]." ";
		if ($b[5] eq "Gamma"){
			$band_labels.=$gamma;
		} else {
			$band_labels.=$b[5]." ";
		}
		if ($b[6] eq "|"){
			$band.=", " ;
		} elsif ($b[6] eq "X"){
			close KPT;
			open KPT, ">tsubo_temp_get_kpath";
			print KPT ("$band \n");
			print KPT ("$band_labels \n");
			close KPT;
			return;
		}
	}
}


sub get_kpath_mass{
#get kpath for mass calculation
#hikisuu: kankaku, SC or Cryst
	&get_kpath_kpoints(@_);
	system ("grep -n '! K' tsubo_temp_get_kpath | sed -e 's/://g' > tsubo_temp1");
#tsubo_temp1 wo yomu
	open TEMP, "tsubo_temp1";
	my (@line, @key, @symbol, @connect);
	while ($a=<TEMP>){
		my @b=&splitall($a);
		push (@line, $b[0]);
		push (@key, $b[6]);
		push (@symbol, $b[7]);
		push (@connect, $b[8]);
	}
	close TEMP;
	system ("rm tsubo_temp1");
	for (my $i=0; $i<=$#line; $i++){
		if ($connect[$i] eq "-"){
			system ("grep -A 7 '${key[$i]} ' tsubo_temp_get_kpath >> tsubo_temp1");
			my $b=$line[$i]+8;
			$a=`sed -n ${b}P tsubo_temp_get_kpath`;
			$a=&rmvcrlf($a);
			system ("echo '${a}! ${symbol[$i]}_to_${symbol[$i+1]} ' >> tsubo_temp1");
			$b=$line[$i+1]-8;
			$a=`sed -n ${b}P tsubo_temp_get_kpath`;
			$a=&rmvcrlf($a);
			system ("echo '${a}! ${symbol[$i+1]}_to_${symbol[$i]} ' >> tsubo_temp1");
			system ("grep -B 7 '${key[$i+1]} ' tsubo_temp_get_kpath | head -7 >> tsubo_temp1");
		} else {	
			system ("grep '${key[$i]} ' tsubo_temp_get_kpath >> tsubo_temp1");
		}
	}
	close TEMP;
	system ("rm tsubo_temp_get_kpath");
}




########
# Miscellaneous subroutines
########


### String

#ataerareta gyou no kaigyou ko-do wo kesu
sub rmvcrlf{
	my $a=$_[0];
	$a=<STDIN> if ($a eq "STDIN");
	$a =~ s/\x0A//g;
	$a =~ s/\x0D//g;
	return ($a);
}

#ataerareta gyou wo hairetu to site kaesu
#saidai 3 youso wo kaesu (tatoeba lattice parameter)
sub split3{
	my $a=&rmvcrlf($_[0]);
	$a=" ".$a;
	my @b=split(/\s+/,$a);
	return (@b[1..3]);
}

#ataerareta gyou wo hairetu to site kaesu
#saidai 4 youso wo kaesu
sub split4{
	my $a=&rmvcrlf($_[0]);
	$a=" ".$a;
	my @b=split(/\s+/,$a);
	shift @b;
	my $c=shift @b;
	my $d=shift @b;
	my $e=shift @b;
	my $f=join(" ",@b);
	return ($c,$d,$e,$f);
}

#ataerareta gyou wo hairetu to site kaesu
#subeteno youso wo kaesu
sub splitall{
	my $a=&rmvcrlf($_[0]);
	$a=" ".$a;
	my @b=split(/\s+/,$a);
	shift @b;
	return (@b);
}

### Array


sub site_in_cell{
#Youso wo 3 ko made 0 kara 1 made no aida ni suru
	my @a=@_;
	for (my $i=0; $i<3; $i++){
		$a[$i]*=1;
		while ($a[$i] < 0){
			$a[$i]++ ;
		}
		$a[$i]-=int($a[$i]);
		$a[$i]=&precise($a[$i],0);
		$a[$i]=0 if (&precise($a[$i],1) == 1);
	}
	return (@a);
}

sub site_is_same{
#futatsu no site ga onaji ka hantei
#hikisuu ha 7 ko:
#tolerance, site 1 no zahyou (3 ko), site 2 no zahyou (3 ko) 
#fractional coordinates de ataery
#onaji nara 1, chigau nara 0 
	my @a=@_;
	my $threshold=shift @a;
	my $b;
	my $c=0;
	for (my $i=0; $i<3; $i++){
		$b=abs ($a[$i]-$a[$i+3]);
		while ($b > 0.8){
			$b--;
		}
		$c+=$b*$b;
	}
	my $ret=0;
	$ret=1 if ($c < $threshold);
	return ($ret);
}

#3x3 gyouretu wo hikisuu ka file kara nyuusyu
sub get_transformation_matrix{
	my @a=@_;
	@a=&args_from_file($_[0]) if ($#a == 0);
	if ($#a == 2){
		return ($a[0],0,0, 0,$a[1],0, 0,0,$a[2], 0,0,0);
	} elsif ($#a == 3){
		return (1,0,0,0,1,0,0,0,1,@a[1..3]);
	} elsif ($#a == 8){
		return (@a, 0,0,0);
	} elsif ($#a == 11){
		return (@a);
	} else {
		die ("Matrix no hikisuu ga okasii desu\n")
	}
}

# get key to sort order
sub order_key_cmp{
	my @orig=@_;
	my @new = sort { $orig[$a] cmp $orig[$b]} 0 .. $#orig;
	return (@new);
}
sub order_key_number{
	my @orig=@_;
	my @new = sort { $orig[$a] <=> $orig[$b]} 0 .. $#orig;
	return (@new);
}


sub expand_array{
#hikisuu wo hairetu wo tenkai
#.. ha jissai no kazu ni tenkai suru
	my @in=@_;
	my @out;
	while (my $a=shift(@in)){
		my @b=split (/\.\./, $a);
		$b[1]=$b[0] if ($b[1] eq "");
		$b[1]=$b[0] if ($b[0]>$b[1]);
		@out = (@out, $b[0] .. $b[1]);
	}
#chouhuku wo kesu
	@out=&unique_array(@out);
#sort
	my @order_key=&order_key_number(@out);
	@out=@out[@order_key];
	return (@out);
}

sub unique_array{
#hairetu de chouhuku suru youso wo kesu
#kekka ha sort sinai
	my @array=@_;
	my %count;
	@array = grep( !$count{$_}++, @array) ;
	return (@array);
}

sub same_array{
#hairetu ga onaji ka douka kakunin (ne wo siyou, mojiretu wo soutei)
#onaji = 0, chigau = 1
	my $n=($#_+1)/2;
	die ("&same_array no hikisuu ga kisuu ko : @_ \n") if (int($n) != $n);
	my $check=0;
	for (my $i=0; $i<$n;$i++){
		$check++ if ($_[$i+$n] ne $_[$i]);
	}
	$check=1 if ($check>1);
	return ($check);
}


#saidaichi
sub max{
	my $m=$_[0];
	for (my $i=1;$i<=$#_;$i++){
		$m=$_[$i] if ($_[$i] > $m);
	}
	return $m;
}

#saishouchi
sub min{
	my $m=$_[0];
	for (my $i=1;$i<=$#_;$i++){
		$m=$_[$i] if ($_[$i] < $m);
	}
	return $m;
}

# within: hikisuu ($a, $b, $c) de
# $b<=$a katu $a<=$c nara 1 soreigai 0

sub within{
	return 0 if ($_[0] < $_[1]);
	return 0 if ($_[0] > $_[2]);
	return 1;
}


sub part_sum_array{
#nyuuryoku ha ($n, @array)
#@array no $n bannme made no atai wo tasu
#tatoeba $n=2 @array=(4,5,6,7,8)
#4+5+6=15 wo kaesu
#konobaai $n no saidaichi ha 4
#nanode $n ga 5 ijyou nara die
#array nagasa check wa okonau
	my @in=@_;
	my $len=shift (@in);
	if (($len < 0) || ($len > $#in)){
		die ("part_sum_array de array nagasa no mondai: $len banme made, hairetu @in\n"); 
	}
	my $sum=0;
	for (my $i=0; $i<=$len; $i++){
		$sum+=$in[$i];
	}
	return ($sum);
}



### Math: vector/matrix operation
##Vector

sub product_vs{
#multiply vector by scalar (scale) 
#scaling factor is last element
#returns array
	my @a;
	for (my $i=0; $i<$#_;$i++){
		$a[$i]=$_[$i]*$_[$#_];
	}
	return @a;
}

sub sum_vv{
#input (@array1, @array2)
#return @array1+@array2
	my $n=($#_+1)/2;
	die ("&sum_vv no hikisuu ga kisuu ko : @_ \n") if (int($n) != $n);
	my @a;
	for (my $i=0; $i<$n;$i++){
		$a[$i]=$_[$i+$n]+$_[$i];
	}
	return (@a);
}

sub diff_vv{
#input (@array1, @array2)
#return @array1-@array2
	my $n=($#_+1)/2;
	die ("&diff_vv no hikisuu ga kisuu ko : @_ \n") if (int($n) != $n);
	my @a;
	for (my $i=0; $i<$n;$i++){
		$a[$i]=$_[$i]-$_[$i+$n];
	}
	return (@a);
}
sub product_vv{
#gets inner product based on two vectors
#returns scalar
	my $n=($#_+1)/2;
	die ("&product_vv no hikisuu ga kisuu ko : @_ \n") if (int($n) != $n);
	my $a;
	for (my $i=0; $i<$n;$i++){
		$a+=$_[$i+$n]*$_[$i];
	}
	return $a;
}

sub cross_vv{
#gets cross product based on two vectors
#returns array
	my $n=($#_+1)/2;
	die ("&cross_vv no hikisuu ga 6 ko de nai : @_ \n") if ($#_ != 5);
	my @a=($_[1]*$_[5]-$_[2]*$_[4],$_[2]*$_[3]-$_[0]*$_[5],$_[0]*$_[4]-$_[1]*$_[3]);
	return @a;
}
sub norm{
#norm of array
	my $a;
	for (my $i=0; $i<=$#_; $i++){
		$a+=$_[$i]*$_[$i];
	}
	$a=sqrt($a);
	return ($a);
}

sub unit_v{
#input @array
#return array with length 1
	my @a=@_;
	my $n=&norm(@a);
	for (my $i=0; $i<=$#a;$i++){
		$a[$i]/=$n;
	}
	return (@a);
}

sub angle_vv{
#angle between (@_[3..5]-@_[0..2]) and (@_[6..8]-@_[0..2])
	my @vec1=&diff_vv(@_[3..5],@_[0..2]);
	my @vec2=&diff_vv(@_[6..8],@_[0..2]);
	my $norm1=&norm(@vec1);
	my $norm2=&norm(@vec2);
	die ("angle_vv: vector 1 ( @_[3..5] ) to origin ( @_[0..2] )  ga onaji") if (&precise($norm1,0,1E-4)==0);
	die ("angle_vv: vector 2 ( @_[6..8] ) to origin ( @_[0..2] ) ga onaji")  if (&precise($norm2,0,1E-4)==0);
	my $angle=&acos(&product_vv(@vec1,@vec2)/$norm1/$norm2);
	return ($angle);
}

##Matrix

sub product_mm{
#multiply matrix
	my @a;
	$a[0]=&product_mm_element(@_,0,1,2,9,12,15);
	$a[1]=&product_mm_element(@_,0,1,2,10,13,16);
	$a[2]=&product_mm_element(@_,0,1,2,11,14,17);
	$a[3]=&product_mm_element(@_,3,4,5,9,12,15);
	$a[4]=&product_mm_element(@_,3,4,5,10,13,16);
	$a[5]=&product_mm_element(@_,3,4,5,11,14,17);
	$a[6]=&product_mm_element(@_,6,7,8,9,12,15);
	$a[7]=&product_mm_element(@_,6,7,8,10,13,16);
	$a[8]=&product_mm_element(@_,6,7,8,11,14,17);
	return @a;
}

sub product_mm_element{
#gets inner product based on two 3x3 matrices 
#first matrix is $_[0] to $[8]
#second matrix is $_[9] to $[17]
#returns the following inner product
	my $i=$_[$_[18]]*$_[$_[21]]+$_[$_[19]]*$_[$_[22]]+$_[$_[20]]*$_[$_[23]];
	return $i;
}

sub invmatrix3{
#gets inverse of 3x3 matrix
	my $d=&det3(@_);
	die ("determinant 0 in @_ \n") if (abs($d) < 1E-12);
	my @a;
	$a[0]=&det2(@_,4,5,7,8)/$d;
	$a[1]=-1*&det2(@_,1,2,7,8)/$d;
	$a[2]=&det2(@_,1,2,4,5)/$d;
	$a[3]=-1*&det2(@_,3,5,6,8)/$d;
	$a[4]=&det2(@_,0,2,6,8)/$d;
	$a[5]=-1*&det2(@_,0,2,3,5)/$d;
	$a[6]=&det2(@_,3,4,6,7)/$d;
	$a[7]=-1*&det2(@_,0,1,6,7)/$d;
	$a[8]=&det2(@_,0,1,3,4)/$d;
	for (my $i=0; $i<=8; $i++){
		$a[$i]=&precise($a[$i],0);
	}
	return @a;
}

sub tmatrix3{
#gets 3x3 transposed matrix
	my @d=($_[0],$_[3],$_[6],$_[1],$_[4],$_[7],$_[2],$_[5],$_[8]);
	return @d;
}

sub det2{
#gets 2x2 determinant from 3x3 matrix
#first 9 items are $_[0] to $[8]
#returns determinant for remaining 4 components
	my $d=$_[$_[9]]*$_[$_[12]]-$_[$_[10]]*$_[$_[11]];
	return $d;
}

#get determinant of 3x3 matrix
sub det3{
	my $d=&product_vv(@_[0..2],&cross_vv(@_[3..8]));
	return $d;
}


sub product_vm{
#gets product of 3-vector and 3x3 matrix
#vector is $_[0] to $[2]
#matrix is $_[3] to $[11]
#returns  array
	my @a;
	$a[0]=$_[0]*$_[3]+$_[1]*$_[6]+$_[2]*$_[9];
	$a[1]=$_[0]*$_[4]+$_[1]*$_[7]+$_[2]*$_[10];
	$a[2]=$_[0]*$_[5]+$_[1]*$_[8]+$_[2]*$_[11];
	return @a;
}

sub product_mv{
#gets product of 3x3 matrix and 3-vector
#vector is $_[9] to $[11]
#matrix is $_[0] to $[8]
#returns array
	my @a;
	$a[0]=$_[0]*$_[9]+$_[1]*$_[10]+$_[2]*$_[11];
	$a[1]=$_[3]*$_[9]+$_[4]*$_[10]+$_[5]*$_[11];
	$a[2]=$_[6]*$_[9]+$_[7]*$_[10]+$_[8]*$_[11];
	return @a;
}

sub apply_symop{
#apply sumop to $_[0] kara $_[2]
#matrix is $_[3] to $[11]
#vector is $_[12] to $[14]
#returns array
	my @a=&product_mv(@_[3..11],@_[0..2]);
	@a=&sum_vv(@a,@_[12..14]);
	return @a;
}


#Rounding

sub floor{
	my $m=$_[0];
	my $tolerance=1E-12;
	if ($m >= 0){
		$m=int($m);
	} else {
		$m=int($m+$tolerance)-1;
	}
	return $m;
}

sub floor_array{
	my @a;
	for (my $i=0;$i<=$#_;$i++){
		$a[$i]=&floor($_[$i]);
	}
	return @a;
}

sub ceiling{
	my $m=$_[0];
	my $tolerance=1E-12;
	if ($m <= 0){
		$m=int($m);
	} else {
		$m=int($m-$tolerance)+1;
	}
	return $m;
}

sub round{
#round to closest integer
	my $m=$_[0];
	if($m < 0){
		$m=$_[0]-0.5;
	} else {
		$m=$_[0]+0.5;
	}
	$m=int($m);
	return $m;
}

#array no subete no youso wo round
sub round_array{
	my @a=@_;
	for (my $i=0;$i<=$#_;$i++){
		$a[$i]=&round($a[$i]);
	}
	return @a;
}

#Trigonometry

sub sind{
	 my $t=sin(atan2(1,1) * $_[0] / 45);
	 return $t;
}

sub cosd{
	 my $t=cos(atan2(1,1) * $_[0] / 45);
	 return $t;
}

sub tand{
	 my $t=&sind($_[0])/&cosd($_[0]);
	 return $t;
}

sub acos{
#arccos
#cos(t) ga wakatteiru baai:
#tan(t) = sqrt(1/cos2(t)-1)*sgn(cos(t))
#cos(t)=0 no baai: acos = 90 
#0-180 degree (0,180 inclusive)
	my $tolerance=1E-12;
	return (90) if (&precise($_[0],0,$tolerance) == 0);
	my $a=1/$_[0]/$_[0]-1;
	$a=0 if (&precise($a,0,$tolerance) == 0);
	$a=sqrt($a);
	if ($_[0] > 0){
		$a=atan2($a,1)/atan2(1,1)/4*180;
	} else {
		$a=atan2($a,-1)/atan2(1,1)/4*180;
	}
	return $a;
}

#Integer/non-array
#seisuu ka douka kakunin
#chigau nara "X" wo kaesu
sub is_integer{
	my $int=&round($_[0]);
	my $err=abs(&precise($_[0]-$int,0,1E-12));
	return "X" if ($err > 1E-12);
	return $int;
}

sub oddeven{
#odd:1 even:0
	my $in=&is_integer($_[0]);
	return 999 if ($in eq "X");
	my $out=($in % 2);
	return $out;
}

#Cleaning
# precise: hikisuu $_[0] ga $_[1] ni chikai baai $_[1] wo kaesu
# tatoeba &precise(1E-17, 0) ha 0 wo kaesu
# $b<=$a katu $a<=$c nara 1 soreigai 0

sub precise{
	my $a=$_[0];
	my $tolerance=1E-12;
	$tolerance=$_[2] if ($_[2] ne "");
	$a=$_[1] if (abs($_[0]-$_[1]) < $tolerance);
	return $a;
}

#integer nara genmitsu ni seisuu ni suru
sub clean_integer{
	my $in=&is_integer($_[0]);
	return ($_[0]) if ($in eq "X");
	return $in;
}

#Integer/array
sub gcmmany{
#ooku no kazu no gcm wo motomeru
#consider negative numbers
#return negative number if all numbers are negative
#return 0 if all numbers are zero
	my @a=@_;
#remove 0
	my @b;
	my $no_positive_flag=1;
	for (my $i=0; $i<=$#a; $i++){
		if ($a[$i] > 0){
			push (@b, $a[$i]);
			$no_positive_flag=0;
		} elsif ($a[$i] < 0) {
			push (@b, -$a[$i]);
		}
	}
#no positive number
	return 0 if ($b[0] eq "");
#one number
	if ($b[1] eq ""){
		$b[0]*=-1 if ($no_positive_flag == 1);
		return ($b[0]) ;
	}
#standard gcm
	my $new=shift(@b);
	while (my $x = shift (@b)){
		$new=&gcmpair($x, $new);
	}
	$new*=-1 if ($no_positive_flag == 1);
	return ($new);
}


sub gcmpair{
#get gcm of pair of numbers
	my $a=&round($_[0]);
	my $b=&round($_[1]);
	die ("gcmpair: $_[0] ga seisuu de nai \n") if (&precise($a-$_[0],0) != 0);
	die ("gcmpair: $_[1] ga seisuu de nai \n") if (&precise($b-$_[1],0) != 0);
	my $flag=0;
	$a=1 if ($a==0);
	$b=1 if ($b==0);
	while(){
		my $c=$a % $b;
		return($b) if ($c == 0);
		$a=$b;
		$b=$c;
	}
}

#array wo GCM de waru
#$_[0] ga "keep_sign" no baai, hu no kazu de waranai
sub minimum_array{
	my $signflag=0;
	if ($_[0] eq "keep_sign"){
		shift @_;
		$signflag=1;
	}
	my @a=@_;
	for (my $i=0;$i<=$#_;$i++){
		die ("minimum array de mondai: $a[$i] ga seisuu de nai\n") if ($a[$i] != (&round($a[$i])));
	}
	my $gcm=&gcmmany(@a);
	$gcm=abs($gcm) if ($signflag==1);
	for (my $i=0;$i<=$#_;$i++){
		$a[$i]/=$gcm;
	}
	return @a;
}




#Misc.

#bunsuu no shousuu, bunsi, bunbo wo kaesu
sub frac2dec{
	my @a=split(/\//,$_[0]);
	die ("frac2dec: $_[0] ga seisuu ka seisuu no bunsuu de nai\n") if ($a[0] !~ /^[-]?[0-9]+$/);
	$a[1]=1 if ($a[1] eq "");
	die ("frac2dec: $_[0] no bunbo ga 0\n") if ($a[1] == 0);
	die ("frac2dec: $_[0] ga seisuu ka seisuu no bunsuu de nai\n") if ($a[1] !~ /^[-]?[0-9]+$/);
	return ($a[0]/$a[1],$a[0],$a[1]);
}

#bunsuu no hairetu wo shyousuu ni nasu
sub frac2dec_array{
	my @a;
	for (my $i=0; $i<=$#_; $i++){
		my @b=&frac2dec($_[$i]);
		$a[$i]=$b[0];
	}
	return (@a);
}

#bunsuu no hairetu wo coprime ni suru
sub frac2dec_gcm{
	my (@denom, @num);
#bunsuu
	for (my $i=0; $i<=$#_; $i++){
		my @b=&frac2dec($_[$i]);
		$num[$i]=$b[1];
		$denom[$i]=$b[2];
	}
#tsubun
	for (my $i=0; $i<=$#_; $i++){
		for (my $j=0; $j<=$#_; $j++){
			$num[$i]*=$denom[$j] if ($i != $j);
		}
	}
	my @a=&minimum_array("keep_sign",@num);
	return (@a);
}

sub projection_vector3D{
#projection of *3D* vector @_[3..5] on @_[0..2]
	my @a=&product_vs(@_[0..2],&product_vv(@_[0..5])/&product_vv(@_[0..2],@_[0..2]));
	return @a;
}

sub projection_normal_vector3D{
#projection of *3D* vector @_[3..5] on PLANE NORMAL TO @_[0..2]
	my @a=&diff_vv(@_[3..5],&projection_vector3D(@_[0..5]));
	return @a;
}

sub product_symop{
#symmetry operation @[0..11] to @[12..23] wo kakezan
#
	my @a=(&product_mm(@_[0..8],@_[12..20]),&sum_vv(&product_mv(@_[0..8],@_[21..23]),@_[9..11]));
	return (@a);
}





#@_ hairetu no aida de saidai no kankaku wo motomeru
#chusin, kankaku wo kaesu
#tatoeba 0, 0.2, 0.8, 1 nara 0.5, 0.6 wo kaesu
#tadashi gap 0.25 katsu center ga 0.625-epsilon ika ga hutatu aru baai ryouhou kaesu

sub max_gap{
	die ("&max_gap de hairetu @_ no nagasa ga 2 miman \n") if ($#_ <= 0);
#unique
	my @array=&unique_array(@_);
#sort
	my @order_key=&order_key_number(@array);
	@array=@array[@order_key];
#	print ("max_gap wo sagasu: @array \n");
	die ("&max_gap no unique hairetu @array no nagasa ga 2 miman \n") if ($#array <= 0);
	my $a=$array[1];
	my $gap=$array[1]-$array[0];
	my $center=($array[1]+$array[0])/2;
	my @return;
#quarter flag: gap ga zenbu 0.25 no baai ha kiken
#alternative no baai center ga 0.375 to 0.625-epsilon wo kaesu!
	if ((abs($array[1]-$array[0]-0.25) < 1E-5) && (abs($array[2]-$array[1]-0.25) < 1E-5)){
		print ("center wo 2ko kaesu: @array \n");
		return (($array[0]+$array[1])/2,0.25,($array[1]+$array[2])/2,0.25) 	}
	for (my $i=2; $i <= $#array; $i++){

		if (($array[$i]-$a) > $gap){
			$gap=$array[$i]-$a;
			$center=($array[$i]+$a)/2;
		}
		$a=$array[$i];
	}
	$center-=0.5 if ($center > 0.5);
	return ($center, $gap);
}



sub gaussian_reduction_2d{
#Gaussian lattice reduction
# http://archive.is/eXKur wo kaihen
#Input: @latvec1,@latvec2,@direction1,@direction2
#Output: reduced
	my @vecu=@_[0..2];
	my @vecv=@_[3..5];
	my @hklu=@_[6..8];
	my @hklv=@_[9..11];
	my @w;
	if (&product_vv(@vecu,@vecv) < 0){
		@vecv=&product_vs(@vecv,-1);
		@hklv=&product_vs(@hklv,-1);
	}
	if (&norm(@vecu) < &norm(@vecv)){
		@w=@vecu;
		@vecu=@vecv;
		@vecv=@w;
		@w=@hklu;
		@hklu=@hklv;
		@hklv=@w;
	}
	while (&norm(@vecu) > &norm(@vecv)){
		my $t=&floor(&product_vv(@vecu,@vecv)/&product_vv(@vecv,@vecv));
		@w=&diff_vv(@vecu,&product_vs(@vecv,$t));
		@vecu=@vecv;
		@vecv=@w;
		@w=&diff_vv(@hklu,&product_vs(@hklv,$t));
		@hklu=@hklv;
		@hklv=@w;
	}
	@w=@vecu;
	@vecu=@vecv;
	@vecv=@w;
	@w=@hklu;
	@hklu=@hklv;
	@hklv=@w;

#kokode 0<&product_vv(@vec1,@vec2)<&norm(@vec1)
#u=@vec2,v=@vec1, v<=u 
#norm(u-v)<v nara kitei torikae
	my (@vec1, @vec2, @hkl1, @hkl2);
	if ((&product_vs(@vecu,@vecv)) > (&norm(@vecv)/2) && (&norm(&diff_vv(@vecu,@vecv))) < &norm(@vecv)){
		@vec1=&diff_vv(@vecu,@vecv);
		@hkl1=&diff_vv(@hklu,@hklv);
		@vec2=&product_vs(@vecv,-1);
		@hkl2=&product_vs(@hklv,-1);
	} else {
		@vec1=@vecv;
		@hkl1=@hklv;
		@vec2=@vecu;
		@hkl2=@hklu;
	}
	return (@vec1,@vec2,@hkl1,@hkl2);
}


#polar: vector r=@_[0..2] no nagasa to, vector a=@_[3..5], b=@_[6..8] no nasu heimen de no @a wo kakudo (hogaku) wo kaesu
#a ga kakudo 0, b ga kakudo 0~180 no aida
#0 kara 360 made no atai wo kaesu
sub polar{
	die ("polar: atai ga 9 yori chiisai: @_ \n") if ($_[8] eq "");
	die ("polar: atai ga 9 yori ookii: @_ \n") if ($_[9] ne "");
	my @r=@_[0..2];
	my @a=@_[3..5];
	my @b=@_[6..8];
# onaji heimen?
	die ("polar: onaji heimen de nai: @_ \n") if ((&det3($_)) != 0);
# nagasa
	my $norm=&norm(@r);
# kakudo
	my $ab=&angle_vv(0,0,0,@a,@b);
	my $ra=&angle_vv(0,0,0,@r,@a);
	my $rb=&angle_vv(0,0,0,@r,@b);
	die ("polar: saigo no vector 2 ga heikou \n") if (&precise($ab,0) == 0);

# shogen 1
	return ($norm,$ra) if (&precise($ra+$rb-$ab,0,1E-5)==0);
# shogen 2
	return ($norm,$ra) if (&precise($ra-$rb-$ab,0,1E-5)==0);
# shogen 3
	return ($norm,$ab+$rb) if (&precise($ra+$rb+$ab,360,1E-5)==360);
# shogen 4
	return ($norm,360-$ra) if (&precise($ra-$rb+$ab,0,1E-5)==0);
	die ("polar: vector @r no vector a @a b @b de kakudo ga hantei dekinai: angle ab = $ab , ra = $ra , rb = $rb \n");
}


#renbunsuu de syousuu wo bunsuu ni suru
sub to_fraction{
	my $num=$_[0];
	my $initial_eps=1E-5;
#integer check
	my $int=&round($num);
	my $dec=&precise($num-$int,0,$initial_eps);
	return ($int) if ($dec == 0);
#negative?
	my $sign=1;
	if ($num < 0){
		$num*=-1;
		$sign=-1;
	}
#renbunsuu tenkai
	my $expansion_eps=1E-3;
	my @expansion;
#seisuu bubun
	$expansion[0]=&floor($num);
	$num-=$expansion[0];
#bunsuu
	my $count=0;
	while ($num >= 0 && $count<10){
		$num=1/$num;
		$count++;
#seisuu ni chikai to seisuu atsukai
		my $round_num=&round($num);
		my $error=&precise($num-&round($num),0,$expansion_eps);
		$num=&round($num) if ($error == 0);
		push (@expansion, int($num));
		last if ($error == 0);
		$num-=int($num);
	}
#@fraction=(bunsi, bunbo)
#saigo kara kumitate
	my @fraction=($expansion[$#expansion],1);

	for (my $i=$#expansion-1; $i>=0;$i--){
		@fraction=($fraction[0]*$expansion[$i]+$fraction[1],$fraction[0]);
	}
	if ($sign == 1){
		return ($fraction[0]."/".$fraction[1]."(".$_[0].")");
	} else {
		return ("-".$fraction[0]."/".$fraction[1]."(".$_[0].")");
	}
}

sub leastsquares{
# https://sci-pursuit.com/math/statistics/least-square-method.html
#input (@array1, @array2)
#return a,b in y=a*x+c
	my $n=($#_+1)/2;
	die ("&sum_vv no hikisuu ga kisuu ko : @_ \n") if (int($n) != $n);
	my @a=@_[0..($n-1)];
	my @b=@_[$n..$#_];
	my ($xave, $yave, $p,$q);
	for (my $i=0; $i<$n; $i++){
		$xave+=$a[$i];
		$yave+=$b[$i];
	}
	$xave/=$n;
	$yave/=$n;
	for (my $i=0; $i<$n; $i++){
		$p+=($a[$i]-$xave)*($b[$i]-$yave);
		$q+=($a[$i]-$xave)*($a[$i]-$xave);
	}
	return ($p/$q, $yave-$p/$q*$xave);
}


#spherical triangle no menseki
#hikisuu ha a,b,c RADIANS!
#S no hanni ha 0<S<2pi thus 0<S/2<pi

sub spherical_triangle_area{
	my $t=(1+$_[0]+$_[1]+$_[2])/sqrt(2*(1+$_[0])*(1+$_[1])*(1+$_[2]));
#S/2=acos($t) https://arxiv.org/pdf/1404.6592.pdf
#S/2=pi/2 (90deg) nara $t=0,S=pi
#use acos(X)=2*atan(sqrt((1-X)/(1+X)))
#thus S=4*atan(sqrt((1-X)/(1+X)))
#$a no hanni ha 0<a
	$t=&precise($t,1);
	my $a=sqrt((1-$t)/(1+$t));
	my $b=4*atan2($a,1);
	return ($b);
}


############
#
#  Lookup table type
#
############

#unique nonpolar orientation = 1, sore igai 0
#input ($space_group_number, $h, $k, $k)
# teigi wa Phys Rev Mater 2, 124603 (2018)
sub is_unique_nonpolar{
	my $sg=$_[0];
#seisuu igai ga areba hajiku
	my $h=$_[1];
	my $k=$_[2];
	my $l=$_[3];
	die ("H $h not integer") if ($h !~ /^[+-]?\d+$/ );
	die ("K $k not integer") if ($k !~ /^[+-]?\d+$/ );
	die ("L $l not integer") if ($l !~ /^[+-]?\d+$/ );
#coprime + all non-negative igai wo hajiku 
	return (0) if (&gcmmany(@_[1..3]) != 1);
#h>0
#	return (0) if ($h < 0);
#individual
	if ($sg == 2){
		return (1) if ($h > 0);
		return (1) if (($h == 0) && ($k >= 0));
	} elsif (&within($sg,3,5)){
		return (1) if (($h >= 0) && ($k == 0));
	} elsif (&within($sg,6,9)){
		return (1) if (! &same_array(@_[1..3],0,1,0));
	} elsif (&within($sg,10,15)){
		return (1) if (($h >= 0) && ($k >= 0));
	} elsif (&within($sg,16,24)){
		return (1) if (($h >= 0) && ($k > 0) && ($l == 0));
		return (1) if (($h > 0) && ($k == 0) && ($l >= 0));
		return (1) if (($h == 0) && ($k >= 0) && ($l > 0));
	} elsif (&within($sg,25,46)){
		return (1) if (($h >= 0) && ($k >= 0) && ($l == 0));
	} elsif (&within($sg,47,74)){
		return (1) if (($h >= 0) && ($k >= 0) && ($l >= 0));
	} elsif (&within($sg,75,80)){
		return (1) if (($h >= 0) && ($k > 0) && ($l == 0));
	} elsif (&within($sg,81,82)){
		return (1) if (($h >= 0) && ($k > 0) && ($l == 0));
		return (1) if (! &same_array(@_[1..3],0,0,1));
	} elsif (&within($sg,83,88)){
		return (1) if (($h >= 0) && ($k > 0) && ($l >= 0));
		return (1) if (! &same_array(@_[1..3],0,0,1));
	} elsif (&within($sg,89,98)){
		return (1) if (($h >= $k) && ($k > 0) && ($l == 0));
		return (1) if (($h > 0) && ($k == 0) && ($l >= 0));
		return (1) if (($h >= 0) && ($k == $h) && ($l > 0));
	} elsif (&within($sg,99,110)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l == 0));
	} elsif (&within($sg,111,114) || &within($sg,121,122)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l == 0));
		return (1) if (($h >= 0) && ($k == 0) && ($l > 0));
	} elsif (&within($sg,115,120)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l == 0));
		return (1) if (($h >= 0) && ($k == $h) && ($l > 0));
	} elsif (&within($sg,123,142)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l >= 0));
	} elsif ($sg == 147){
#2019/11/13 ERROR FOUND IN PRM 2 124603
#ORIGINAL
#		return (1) if (($h >= 0) && ($k >= 0));
#NEW
#h>0,k>=0,l
#0,0,1
		if (($h == 0) && ($l == 0)){
			return (1) if ($l == 1);
			return (0);
		}
		return (1) if (($h > 0) && ($k >= 0));
	} elsif ($sg == 148){
#2019/11/13 ERROR FOUND IN PRM 2 124603
#ORIGINAL
#		my @r=&product_mv(2,1,1,-1,1,1,-1,-2,1,@_[1..3]);
#		return (1) if (($r[0] > $r[2]) && ($r[1] > $r[2]) && ($r[0] > 0) && ($r[1] > 0) && ($r[0] != $r[1]));
#		return (1) if (($r[0] > 0) && ($r[1] == $r[0]) && ($r[2] != 0));
#		return (1) if (($r[0] > 0) && ($r[2] == 0));
#NEW
#FIRST 2 OK
#h>=0,k!=0,l=0
#1,0,0
		my @r=&product_mv(2,1,1,-1,1,1,-1,-2,1,@_[1..3]);
		return (1) if (($r[0] > $r[2]) && ($r[1] > $r[2]) && ($r[0] > 0) && ($r[1] > 0) && ($r[0] != $r[1]));
		return (1) if (($r[0] > 0) && ($r[1] == $r[0]) && ($r[2] != 0));
		return (1) if (($r[0] > 0) && ($r[1] != 0) &&($r[2] == 0));
		return (1) if (($r[0] == 1) && ($r[1] == 0) &&($r[2] == 0));
	} elsif (($sg == 149) || ($sg == 151) || ($sg == 153)){
		return (1) if (($h >= 0) && ($k == $h));
	} elsif ($sg == 155){
		my @r=&product_mv(2,1,1,-1,1,1,-1,-2,1,@_[1..3]);
		return (1) if (($r[0] >= 0) && ($r[1] == $r[0]));
	} elsif (($sg == 150) || ($sg == 152) || ($sg == 154)){
		return (1) if (($h >= 0) && ($k == 0));
	} elsif (($sg == 156) || ($sg == 158)){
		return (1) if (! &same_array(@_[1..3],1,1,0));
	} elsif (($sg == 157) || ($sg == 159)){
		return (1) if (! &same_array(@_[1..3],1,0,0));
	} elsif (&within($sg,160,161)){
		my @r=&product_mv(2,1,1,-1,1,1,-1,-2,1,@_[1..3]);
		return (1) if (! &same_array(@r,1,0,-1));
	} elsif (&within($sg,162,163)){
		return (1) if (($h >= $k) && ($k > 0));
		return (1) if (($h >= 0) && ($k == 0) && ($l >= 0));
	} elsif (&within($sg,164,165)){
#2019/11/13 TYPO FOUND IN PRM 2 124603 should be -3m1
		return (1) if (($h > $k) && ($k >= 0));
		return (1) if (($h >= 0) && ($k == $h) && ($l >= 0));
	} elsif (&within($sg,166,167)){
		my @r=&product_mv(2,1,1,-1,1,1,-1,-2,1,@_[1..3]);
		return (1) if (($r[0] >= $r[1]) && ($r[1] >= $r[2]) && ($r[1] > 0) && ($r[2] != 0));
		return (1) if (($r[0] >= abs($r[1])) && (abs($r[1]) >= 0) && ($r[2] == 0));
	} elsif (&within($sg,168,173)){
		return (1) if (($h > 0) && ($k >= 0) && ($l == 0));
	} elsif ($sg == 174){
		return (1) if (! &same_array(@_[1..3],0,0,1));
	} elsif (&within($sg,175,176)){
		return (1) if (($h > 0) && ($k >= 0) && ($l >= 0));
		return (1) if (! &same_array(@_[1..3],0,0,1));
	} elsif (&within($sg,177,182)){
		return (1) if (($h >= $k) && ($k > 0) && ($l == 0));
		return (1) if (($h > 0) && ($k == 0) && ($l >= 0));
		return (1) if (($h >= 0) && ($k == $h) && ($l > 0));
	} elsif (&within($sg,183,186)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l == 0));
	} elsif (&within($sg,187,188)){
		return (1) if (($h >= 0) && ($k == $h) && ($l >= 0));
	} elsif (&within($sg,189,190)){
		return (1) if (($h >= 0) && ($k == 0) && ($l >= 0));
	} elsif (&within($sg,190,194)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l >= 0));
	} elsif (&within($sg,195,199)){
		return (1) if (($h >= 0) && ($k >= 0) && ($l -= 0));
	} elsif (&within($sg,200,206)){
		return (1) if (($h > $l) && ($k > $l) && ($l >= 0) && ($h != $k));
		return (1) if (($h >= 0) && ($k == $h) && ($l >= 0));
	} elsif (&within($sg,207,214)){
		return (1) if (($h >= 0) && ($k == $h) && ($l >= 0));
		return (1) if (($h > $k) && ($k > 0) && ($l == 0));
	} elsif (&within($sg,215,220)){
		return (1) if (($h >= $k) && ($k >= 0) && ($l == 0));
	} elsif (&within($sg,221,230)){
		return (1) if (($h >= $k) && ($k >= $l) && ($l >= 0));
	}
	return (0);
}



# orientation @_[1..3] ni taisi, unique nonpolar orientation wo kaesu
# input ($space_group_number, $h, $k, $k)
# unique nonpolar orientation no teigi wa Phys Rev Mater 2, 124603 (2018)
sub to_unique_nonpolar{
	my ($sg, $h, $k, $l)=@_[0..3];
	my @ret;
#seisuu igai ga areba hajiku
	die ("H $h not integer") if ($h !~ /^[+-]?\d+$/ );
	die ("K $k not integer") if ($k !~ /^[+-]?\d+$/ );
	die ("L $l not integer") if ($l !~ /^[+-]?\d+$/ );
#coprime to 0 0 0 igai wo hajiku 
	my $flag=&gcmmany($h,$k,$l);
	return ("All_zero") if ($flag == 0);
#coprime + all non-negative igai wo hajiku 
	return ("Not_coprime") if (abs($flag) >= 2);
	($h,$k,$l)=(-$h,-$k,-$l) if ($flag == -1);
#individual
	if ($sg == 2){
		return ($h,$k,$l) if ($h > 0);
		return (-$h,-$k,-$l) if ($h < 0);
		return (0,$k,$l) if ($k > 0);
		return (0,-$k,-$l) if ($k < 0);
		return (0,0,1);
	} elsif (&within($sg,3,5)){
		return ("Nonpolar") if ($k != 0);
		return ($h,0,$l) if ($h > 0);
		return (-$h,0,-$l) if ($h < 0);
		return (0,0,1);
	} elsif (&within($sg,6,9)){
		return (0,1,0) if (! &same_array($h,$l,0,0));
	} elsif (&within($sg,10,15)){
#two zero
		return (1,0,0) if (! &same_array($k,$l,0,0));
		return (0,1,0) if (! &same_array($h,$l,0,0));
		return (0,0,1) if (! &same_array($h,$k,0,0));
		($h,$k,$l)=(-$h,-$k,-$l) if ($h < 0);
		return ($h,abs($k),$l);
	} elsif (&within($sg,16,24)){
		return ("Nonpolar") if (($h*$k*$l) != 0);
		return (abs($h),abs($k),abs($l));
	} elsif (&within($sg,25,46)){
		return ("Nonpolar") if ($l != 0);
		return (abs($h),abs($k),0);
	} elsif (&within($sg,47,74)){
		return (abs($h),abs($k),abs($l));
	} elsif (&within($sg,75,80)){
		return ("Nonpolar") if ($l != 0);
		return ($h,$k,0) if (($h >= 0) && ($k > 0));
		return ($k,-$h,0) if (($h < 0) && ($k >= 0));
		return (-$h,-$k,0) if (($h <= 0) && ($k < 0));
		return (-$k,$h,0) if (($h > 0) && ($k <= 0));
	} elsif (&within($sg,81,82)){
		return (0,0,1) if (! &same_array($h,$k,0,0));
		return ("Nonpolar") if ($l != 0);
		return ($h,$k,0) if (($h >= 0) && ($k > 0));
		return ($k,-$h,0) if (($h < 0) && ($k >= 0));
		return (-$h,-$k,0) if (($h <= 0) && ($k < 0));
		return (-$k,$h,0) if (($h > 0) && ($k <= 0));
	} elsif (&within($sg,83,88)){
		return (0,0,1) if (! &same_array($h,$k,0,0));
		($h,$k,$l)=(-$h,-$k,-$l) if ($l < 0);
		return ($h,$k,$l) if (($h >= 0) && ($k > 0));
		return ($k,-$h,$l) if (($h < 0) && ($k >= 0));
		return (-$h,-$k,$l) if (($h <= 0) && ($k < 0));
		return (-$k,$h,$l) if (($h > 0) && ($k <= 0));
	} elsif (&within($sg,89,98)){
		if ($l == 0){
			if (($h < 0) && ($k >= 0)){
				($h,$k)=($k,-$h);
			} elsif (($h <= 0) && ($k < 0)){
				($h,$k)=(-$h,-$k);
			} elsif (($h > 0) && ($k <= 0)){
				($h,$k)=(-$k,$h);
			}
			($h,$k)=($k,$h) if ($h < $k);
			return ($h,$k,0);
		}
		return (abs($k),0,abs($l)) if ($h == 0);
		return (abs($h),0,abs($l)) if ($k == 0);
		return (abs($h),abs($k),abs($l)) if ((abs($h)) == (abs($k)));
		return ("Nonpolar");
	} elsif (&within($sg,99,110)){
		return ("Nonpolar") if ($l != 0);
		($h,$k)=(abs($h),abs($k));
		($h,$k)=($k,$h) if ($h < $k);
		return ($h,$k,0);
	} elsif (&within($sg,111,114) || &within($sg,121,122)){
		return ("Nonpolar") if (($h*$k*$l) != 0);
		if ($l == 0){
			($h,$k)=(abs($h),abs($k));
			($h,$k)=($k,$h) if ($h < $k);
			return ($h,$k,0);
		} elsif ($h == 0){
			return (abs($k),0,abs($l));
		} else {
			return (abs($h),0,abs($l));
		}
	} elsif (&within($sg,115,120)){
		if ($l == 0){
			($h,$k)=(abs($h),abs($k));
			($h,$k)=($k,$h) if ($h < $k);
			return ($h,$k,0);
		} elsif ((abs($h)) == (abs($k))){
			return (abs($h),abs($h),abs($l));
		} else {
			return ("Nonpolar") ;
		}
	} elsif (&within($sg,123,142)){
		($h,$k,$l)=(abs($h),abs($k),abs($l));
		($h,$k)=($k,$h) if ($h < $k);
		return ($h,$k,$l);
	} elsif ($sg == 147){
#-3 hP: h,k,i no hutatu wo sei ni suru
#h,k,i no zero ha 0,2,3
#zero=3
		return (0,0,1) if (! &same_array($h,$k,0,0));
		my $i=0-$h-$k;
#zero=0: ihutatu sei ni site $i wo fu ni suru
#h,k,i no futatu ga sei
		if (($h*$k) != 0) {
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,-$l);
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,-$l);
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,-$l);
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,$l);
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,$l);
			}
			return ($h,$k,$l);
		}
#zero=1: $h ka $k ga non-zero
		if ($h < 0) {
			($h,$k,$l)=(-$h,0,-$l);
		} elsif ($k > 0){
			($h,$k,$l)=($k,0,$l);
		} elsif ($k < 0){
			($h,$k,$l)=(-$k,0,-$l);
		}
		return ($h,$k,$l);

	} elsif ($sg == 148){
#-3 hR: h,k,i no hutatu wo sei ni suru
		($h,$k,$l)=&product_mv(2,1,1,-1,1,1,-1,-2,1,$h,$k,$l);
#00x,0x0,x00
		return (0,1,0) if (! &same_array($k,$l,0,0));
		return (0,1,0) if (! &same_array($h,$l,0,0));
		return (0,1,0) if (! &same_array($h,$k,0,0));

#nonzero=3
		if (($h*$k*$l) != 0){
#negative 2,3 nara hanten
			if (($h < 0) && ($k < 0) && ($l < 0)){
				($h,$k,$l)=(-$h,-$k,-$l);
			} elsif (($h > 0) && ($k < 0) && ($l < 0)) {
				($h,$k,$l)=(-$h,-$k,-$l);
			} elsif (($h < 0) && ($k > 0) && ($l < 0)) {
				($h,$k,$l)=(-$h,-$k,-$l);
			} elsif (($h < 0) && ($k < 0) && ($l > 0)) {
				($h,$k,$l)=(-$h,-$k,-$l);
			}
#min wo l ni suru
			if (&min($h,$k,$l) == $h) {
				($h,$k,$l)=($k,$l,$h);
			} elsif (&min($h,$k,$l) == $k) {
				($h,$k,$l)=($l,$h,$k);
			}
#hhl, hhh OK
#00l nado ha nai
#hh ga negative nara flip
		} elsif ($h == $k) {
			($h,$k,$l)=(-$h,-$h,-$l) if ($h < 0);
#hkk 
#flip, shuffle
		} elsif ($k == $l){
			if ($k < 0){
				($h,$k,$l)=(-$k,-$k,-$h);
			} else {
				($h,$k,$l)=($k,$k,$h);
			}
#hkh 
		} elsif ($l == $h) {
			if ($k < 0){
				($h,$k,$l)=(-$l,-$l,-$k);
			} else {
				($h,$k,$l)=($l,$l,$k);
			}
		} 
		return ($h,$k,$l);
	} elsif (($sg == 149) || ($sg == 151) || ($sg == 153)){
		my $i=0-$h-$k;
		if ($h == $k) {
			($h,$k,$i,$l)=($h,$h,$i,$l);
		} elsif ($h == $i) {
			($h,$k,$i,$l)=($h,$h,$k,$l);
		} elsif ($k == $i) {
			($h,$k,$i,$l)=($k,$k,$h,$l);
		} else {
			return ("Nonpolar");
		}
		($h,$k,$l)=(-$h,-$k,-$l) if ($h < 0);
		return ($h,$k,$l);
	} elsif ($sg == 155){
		($h,$k,$l)=&product_mv(2,1,1,-1,1,1,-1,-2,1,$h,$k,$l);
		if ($h == $k) {
			($h,$k,$l)=($h,$h,$l);
			($h,$k,$l)=(-$h,-$h,-$l) if ($h < 0);
			return ($h,$k,$l);
		} elsif ($h == $l) {
			($h,$k,$l)=($h,$h,$k);
			($h,$k,$l)=(-$h,-$h,-$l) if ($h < 0);
			return ($h,$k,$l);
		} elsif ($k == $l) {
			($h,$k,$l)=($k,$k,$h);
			($h,$k,$l)=(-$h,-$h,-$l) if ($h < 0);
			return ($h,$k,$l);
		} else {
			return ("Nonpolar");
		}
	} elsif (($sg == 150) || ($sg == 152) || ($sg == 154)){
#00x
		return (0,0,1) if (! &same_array($h,$k,0,0));
#other
		my $i=0-$h-$k;
		if ($h == 0) {
			($h,$k,$i,$l)=($i,0,$k,$l);
			($h,$k,$l)=(-$h,-$k,-$l) if ($h < 0);
			return ($h,$k,$l);
		} elsif ($k == 0) {
			($h,$k,$l)=(-$h,-$k,-$l) if ($h < 0);
			return ($h,$k,$l);
		} elsif ($i == 0) {
			($h,$k,$i,$l)=($k,0,$h,$l);
			($h,$k,$l)=(-$h,-$k,-$l) if ($h < 0);
			return ($h,$k,$l);
		}
	} elsif (($sg == 156) || ($sg == 158)){
		return (1,1,0) if (! &same_array($h,$k,$l,1,1,0));
		return (1,1,0) if (! &same_array($h,$k,$l,-2,1,0));
		return (1,1,0) if (! &same_array($h,$k,$l,1,-2,0));
		return (1,1,0) if (! &same_array($h,$k,$l,-1,-1,0));
		return (1,1,0) if (! &same_array($h,$k,$l,-1,2,0));
		return (1,1,0) if (! &same_array($h,$k,$l,2,-1,0));
	} elsif (($sg == 157) || ($sg == 159)){
		return (1,0,0) if (! &same_array($h,$k,$l,1,0,0));
		return (1,0,0) if (! &same_array($h,$k,$l,-1,1,0));
		return (1,0,0) if (! &same_array($h,$k,$l,0,-1,0));
		return (1,0,0) if (! &same_array($h,$k,$l,0,1,0));
		return (1,0,0) if (! &same_array($h,$k,$l,1,-1,0));
		return (1,0,0) if (! &same_array($h,$k,$l,-1,0,0));
	} elsif (&within($sg,160,161)){
		($h,$k,$l)=&product_mv(2,1,1,-1,1,1,-1,-2,1,$h,$k,$l);
		return (1,0,1) if (! &same_array($h,$k,$l,0,1,-1));
		return (1,0,1) if (! &same_array($h,$k,$l,-1,0,1));
		return (1,0,1) if (! &same_array($h,$k,$l,1,-1,0));
		return (1,0,1) if (! &same_array($h,$k,$l,1,0,-1));
		return (1,0,1) if (! &same_array($h,$k,$l,0,-1,1));
		return (1,0,1) if (! &same_array($h,$k,$l,-1,1,0));

	} elsif (&within($sg,162,163)){
#hki zero=3
		return (0,0,1) if (! &same_array($h,$k,0,0));
		my $i=0-$h-$k;
#zero=0: h,k wo sei ni shite $i wo fu ni suru
#h,k,i no futatu ga sei
		if (($h*$k) != 0) {
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,-$l);
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,-$l);
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,-$l);
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,$l);
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,$l);
			}
			($h,$k,$l)=($k,$h,-$l) if ($h < $k);
		} else {
#zero = 1
			($h,$k,$i,$l)=($i,0,$k,$l) if ($h == 0);
			($h,$k,$i,$l)=($k,0,$h,$l) if ($i == 0);
			($h,$k,$l)=(-$h,0,-$l) if ($h < 0);
		}
		return ($h,$k,$l);
	} elsif (&within($sg,164,165)){
		my $i=0-$h-$k;
#hhl wo tsubusu, 001 mo makizoe
		return (abs($h), abs($h), abs($l)) if ($h == $k);
		return (abs($h), abs($h), abs($l)) if ($h == $i);
		return (abs($k), abs($k), abs($l)) if ($k == $i);
#h0l
#hki no doreka ga 0
		if (($h*$k*$i) == 0){
			($h,$k,$i,$l)=($i,0,$k,$l) if ($h == 0);
			($h,$k,$i,$l)=($k,0,$h,$l) if ($i == 0);
			($h,$k,$l)=(-$h,0,-$l) if ($h < 0);
			return ($h,$k,$l);
		}
#hki zenbu non-zero
		if (($h < 0) && ($k < 0)){
			($h,$k,$i,$l)=(-$h,-$k,-$i,-$l);
		} elsif (($h < 0) && ($i < 0)){
			($h,$k,$i,$l)=(-$i,-$h,-$k,-$l);
		} elsif (($k < 0) && ($i < 0)){
			($h,$k,$i,$l)=(-$k,-$i,-$h,-$l);
		} elsif (($h > 0) && ($i > 0)){
			($h,$k,$i,$l)=($i,$h,$k,$l);
		} elsif (($k > 0) && ($i > 0)){
			($h,$k,$i,$l)=($k,$i,$h,$l);
		}
		($h,$k,$l)=($k,$h,$l) if ($h < $k);
		return ($h,$k,$l);
	} elsif (&within($sg,166,167)){
		($h,$k,$l)=&product_mv(2,1,1,-1,1,1,-1,-2,1,$h,$k,$l);
#no zero
		if (($h*$k*$l) != 0){
#sort, h largest
			($h,$k)=($k,$h) if ($k>$h);
			($k,$l)=($l,$k) if ($l>$k);
			($h,$k)=($k,$h) if ($k>$h);
#mannnaka ga + ni naru youni suru
			($h,$k,$l)=(-$l,-$k,-$h) if ($k < 0);
			return ($h,$k,$l);
		}
#two zero
		return (1,0,0) if (! &same_array($k,$l,0,0));
		return (1,0,0) if (! &same_array($h,$l,0,0));
		return (1,0,0) if (! &same_array($h,$k,0,0));
#one zero
		if ($h == 0) {
			($h,$k,$l)=($k,$l,0);
		} elsif ($k == 0) {
			($h,$k,$l)=($l,$h,0);
		}
		($h,$k,$l)=($k,$h,0) if (abs($k) > abs($h));
		($h,$k,$l)=(-$h,-$k,0) if ($h < 0);
		return ($h,$k,$l);
	} elsif (&within($sg,168,173)){
		return ("Nonpolar") if ($l != 0);
		my $i=0-$h-$k;
#h,k nonzero
		if (($h*$k) != 0){
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,0);
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,0);
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,0);
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,0);
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,0);
			}
			return ($h,$k,$l);
		} 
#one zero
		return (1,0,0);
	} elsif ($sg == 174){
		return (0,0,1) if (! &same_array($h,$k,0,0));
	} elsif (&within($sg,175,176)){
#001
		return (0,0,1) if (! &same_array($h,$k,0,0));
		my $i=0-$h-$k;
#h,k nonzero
		if (($h*$k) != 0){
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,abs($l));
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,abs($l));
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,abs($l));
			} elsif (($h > 0) && ($k > 0)){
				($h,$k,$i,$l)=($h,$k,$i,abs($l));
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,abs($l));
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,abs($l));
			}
			return ($h,$k,$l);
		} 
#one zero
		if ($h == 0) {
			($h,$k,$l)=(0,abs($k),abs($l));
		} elsif ($k == 0) {
			($h,$k,$l)=(0,abs($h),abs($l));
		}
		return ($h,$k,$l);
	} elsif (&within($sg,177,182)){
		my $i=0-$h-$k;
#h,k,l nonzero
		if (($h*$k*$l) == 0){
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,0);
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,0);
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,0);
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,0);
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,0);
			}
			($h,$k,$l)=($k,$h,$l) if ($h < $k);
			return ($h,$k,$l);
		}
#0 in h,k,i, 001 mo makizoe
		if ($h == 0) {
			return (abs($k),0,abs($l));
		} elsif ($k == 0) {
			return (abs($h),0,abs($l));
		}
#hhl
		if ($h == $k) {
			return (abs($h),abs($h),abs($l));
		} elsif ($k == $i) {
			return (abs($k),abs($k),abs($l));
		} elsif ($h == $i) {
			return (abs($h),abs($h),abs($l));
		}
	} elsif (&within($sg,183,186)){
#l=0
		if ($l == 0){
			my $i=0-$h-$k;
#one 0
			if ($h == 0) {
				return (abs($k),0,abs($l));
			} elsif ($k == 0) {
				return (abs($h),0,abs($l));
			}
#no 0
			if (($h < 0) && ($k < 0)){
				($h,$k,$i,$l)=(-$h,-$k,-$i,0);
			} elsif (($h < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$i,-$h,-$k,0);
			} elsif (($k < 0) && ($i < 0)){
				($h,$k,$i,$l)=(-$k,-$i,-$h,0);
			} elsif (($h > 0) && ($i > 0)){
				($h,$k,$i,$l)=($i,$h,$k,0);
			} elsif (($k > 0) && ($i > 0)){
				($h,$k,$i,$l)=($k,$i,$h,0);
			}
			($h,$k,$l)=($k,$h,$l) if ($h < $k);
			return ($h,$k,$l);
		}
	} elsif (&within($sg,187,188)){
		my $i=0-$h-$k;
		if ($h == $k) {
			return (abs($k),abs($k),abs($l));
		} elsif ($k == $i) {
			return (abs($k),abs($k),abs($l));
		} elsif ($h == $i) {
			return (abs($h),abs($h),abs($l));
		}
	} elsif (&within($sg,189,190)){
		my $i=0-$h-$k;
		return (abs($k), 0, abs($l)) if ($h == 0);
		return (abs($h), 0, abs($l)) if ($k == 0);
		return (abs($h), 0, abs($l)) if ($i == 0);
	} elsif (&within($sg,190,194)){
		my $i=0-$h-$k;
#0 in h,k,i
		return (abs($k), 0, abs($l)) if ($h == 0);
		return (abs($h), 0, abs($l)) if ($k == 0);
		return (abs($h), 0, abs($l)) if ($i == 0);
#h,k nonzero
		if (($h < 0) && ($k < 0)){
			($h,$k,$i,$l)=(-$h,-$k,-$i,abs($l));
		} elsif (($h < 0) && ($i < 0)){
			($h,$k,$i,$l)=(-$i,-$h,-$k,abs($l));
		} elsif (($k < 0) && ($i < 0)){
			($h,$k,$i,$l)=(-$k,-$i,-$h,abs($l));
		} elsif (($h > 0) && ($k > 0)){
			($h,$k,$i,$l)=($h,$k,$i,abs($l));
		} elsif (($h > 0) && ($i > 0)){
			($h,$k,$i,$l)=($i,$h,$k,abs($l));
		} elsif (($k > 0) && ($i > 0)){
			($h,$k,$i,$l)=($k,$i,$h,abs($l));
		}
		($h,$k,$l)=($k,$h,$l) if ($h < $k);
		return ($h,$k,$l);
	} elsif (&within($sg,195,199)){
		return ("Nonpolar") if (($h*$k*$l) != 0);
		($h,$k,$l)=(abs($h),abs($k),abs($l));
#sort, h largest
		($h,$k)=($k,$h) if ($k>$h);
		($k,$l)=($l,$k) if ($l>$k);
		($h,$k)=($k,$h) if ($k>$h);
		return ($h,$k,$l);
	} elsif (&within($sg,200,206)){
		($h,$k,$l)=(abs($h),abs($k),abs($l));
#hhl
		return ($h,$h,$l) if ($h == $k);
		return ($h,$h,$k) if ($h == $l);
		return ($k,$k,$h) if ($k == $l);
#hkl subete chigau
		return ($h,$k,$l) if (&min($h,$k,$l) == $l);
		return ($l,$h,$k) if (&min($h,$k,$l) == $k);
		return ($k,$l,$h) if (&min($h,$k,$l) == $h);
	} elsif (&within($sg,207,214)){
		($h,$k,$l)=(abs($h),abs($k),abs($l));
#hhl, 00l mo makizoe
		return ($h,$h,$l) if ($h == $k);
		return ($h,$h,$k) if ($h == $l);
		return ($k,$k,$h) if ($k == $l);
#no hkl
		return ("Nonpolar") if (($h*$k*$l) != 0);
#sort, h largest
		($h,$k)=($k,$h) if ($k>$h);
		($k,$l)=($l,$k) if ($l>$k);
		($h,$k)=($k,$h) if ($k>$h);
		return ($h,$k,$l);
	} elsif (&within($sg,215,220)){
#no hkl
		return ("Nonpolar") if (($h*$k*$l) != 0);
		($h,$k,$l)=(abs($h),abs($k),abs($l));
#sort, h largest
		($h,$k)=($k,$h) if ($k>$h);
		($k,$l)=($l,$k) if ($l>$k);
		($h,$k)=($k,$h) if ($k>$h);
		return ($h,$k,$l);
	} elsif (&within($sg,221,230)){
		($h,$k,$l)=(abs($h),abs($k),abs($l));
#sort, h largest
		($h,$k)=($k,$h) if ($k>$h);
		($k,$l)=($l,$k) if ($l>$k);
		($h,$k)=($k,$h) if ($k>$h);
		return ($h,$k,$l);
	}
	return ("Nonpolar");
}


#Kessyogaku  ni motozuita
#Brillouin zone no ten
#$_[0] ni inversion wo kyousei suru ka douka ireru 
# ("Y"=force inversion symmetry)
#Shuturyoku ha array ($latticetype, @zahyou, $symbol ...)
#@zahyou: eg (0,0,0), $symbol: eg Gamma, A,
#Kyousei teki ni crystallographic primitive ni suru
#tadashi triclinic ha KOKODE reduced ni suru
sub bz_point{
	my $extended_bravais=&get_crystallographic_cell ("conventional");
	my $inversion=substr($extended_bravais,3,3);
	my $bravais=substr($extended_bravais,0,2);
	my $centering=substr($extended_bravais,1,1);
	my $latticetype=substr($extended_bravais,0,3);
	$inversion="Y" if ($_[0] eq "Y");
	my @abcreal=&get_abc(@latvec);
	my ($a, $b, $c)=@abcreal[0..2];
	my $cal=&cosd($abcreal[3]);
	my $cbe=&cosd($abcreal[4]);
	my $cga=&cosd($abcreal[5]);
	my $sal=&sind($abcreal[3]);
	my $sbe=&sind($abcreal[4]);
	my $sga=&sind($abcreal[5]);
	my @t=($extended_bravais,0,0,0,"Gamma");
	if ($bravais eq "cP"){
		push @t, (1/2,1/2,1/2,"R");
		push @t, (1/2,1/2,0,"M");
		push @t, (0,1/2,0,"X");
		push @t, (1/2,0,0,"X1") ;
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,-1/2,"R'");
			push @t, (-1/2,-1/2,0,"M'");
			push @t, (0,-1/2,0,"X'");
			push @t, (-1/2,0,0,"X1'");
		}
	} elsif ($bravais eq "cF"){
		push @t, (1/2,0,1/2,"X");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (1/2,1/4,3/4,"W");
		push @t, (3/4,1/4,1/2,"W2");
		push @t, (3/8,3/8,3/4,"K");
		push @t, (5/8,1/4,5/8,"U");
		if ($inversion ne "Y"){
			push @t, (-1/2,0,-1/2,"X'");
			push @t, (-1/2,-1/2,-1/2,"L'");
			push @t, (-1/2,-1/4,-3/4,"W'");
			push @t, (-3/4,-1/4,-1/2,"W2'");
			push @t, (-3/8,-3/8,-3/4,"K'");
			push @t, (-5/8,-1/4,-5/8,"U'");
		}
	} elsif ($bravais eq "cI"){
		push @t, (1/2,-1/2,1/2,"H");
		push @t, (1/4,1/4,1/4,"P");
		push @t, (0,0,1/2,"N");
		if ($inversion ne "Y"){
			push @t, (-1/2,1/2,-1/2,"H'");
			push @t, (-1/4,-1/4,-1/4,"P'");
			push @t, (0,0,-1/2,"N'");
		}
	} elsif ($bravais eq "tP"){
		push @t, (0,0,1/2,"Z");
		push @t, (1/2,1/2,0,"M");
		push @t, (1/2,1/2,1/2,"A");
		push @t, (0,1/2,1/2,"R");
		push @t, (0,1/2,0,"X");
		if ($inversion ne "Y"){
			push @t, (0,0,-1/2,"Z");
			push @t, (-1/2,-1/2,0,"M");
			push @t, (-1/2,-1/2,-1/2,"A");
			push @t, (0,-1/2,-1/2,"R");
			push @t, (0,-1/2,0,"X");
		}
	} elsif ($latticetype eq "tI1"){
# tI1:c<a
		my $y=(1+$c*$c/$a/$a)/4;
		push @t, (-1/2,1/2,1/2,"M");
		push @t, (0,0,1/2,"X");
		push @t, (1/4,1/4,1/4,"P");
		push @t, ($y,$y,-$y,"Z");
		push @t, (-$y,1-$y,$y,"Z0");
		push @t, (0,1/2,0,"N");
		if ($inversion ne "Y"){
			push @t, (1/2,-1/2,-1/2,"M'");
			push @t, (0,0,-1/2,"X'");
			push @t, (-1/4,-1/4,-1/4,"P'");
			push @t, (-$y,-$y,$y,"Z'");
			push @t, ($y,-1+$y,-$y,"Z0'");
			push @t, (0,-1/2,0,"N'");
		}
	} elsif ($latticetype eq "tI2"){
# tI2:c>a
		my $y=(1+$a*$a/$c/$c)/4;
		my $z=$a*$a/2/$c/$c;
		push @t, (1/2,1/2,-1/2,"M");
		push @t, (0,0,1/2,"X");
		push @t, (1/4,1/4,1/4,"P");
		push @t, (0,1/2,0,"N");
		push @t, (-$y,$y,$y,"S0");
		push @t, ($y,1-$y,-$y,"S");
		push @t, (-$z,$z,1/2,"R");
		push @t, (1/2,1/2,-$z,"G");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,1/2,"M'");
			push @t, (0,0,-1/2,"X'");
			push @t, (-1/4,-1/4,-1/4,"P'");
			push @t, (0,-1/2,0,"N'");
			push @t, ($y,-$y,-$y,"S0'");
			push @t, (-$y,-1+$y,$y,"S'");
			push @t, ($z,-$z,-1/2,"R'");
			push @t, (-1/2,-1/2,$z,"G'");
		}
	} elsif ($bravais eq "oP"){
		push @t, (1/2,0,0,"X");
		push @t, (0,0,1/2,"Z");	
		push @t, (1/2,0,1/2,"U");
		push @t, (0,1/2,0,"Y");
		push @t, (1/2,1/2,0,"S");
		push @t, (0,1/2,1/2,"T");
		push @t, (1/2,1/2,1/2,"R");
		if ($inversion ne "Y"){
			push @t, (-1/2,0,0,"X'");
			push @t, (0,0,-1/2,"Z'");	
			push @t, (-1/2,0,-1/2,"U'");
			push @t, (0,-1/2,0,"Y'");
			push @t, (-1/2,-1/2,0,"S'");
			push @t, (0,-1/2,-1/2,"T'");
			push @t, (-1/2,-1/2,-1/2,"R'");
		}
	} elsif ($latticetype eq "oF1"){
#oF1:a^-2>b^-2+c^-2
		my $z=(1+$a*$a/$b/$b-$a*$a/$c/$c)/4;
		my $y=(1+$a*$a/$b/$b+$a*$a/$c/$c)/4;
		push @t, (1,1/2,1/2,"T");
		push @t, (1/2,1/2,0,"Z");
		push @t, (1/2,0,1/2,"Y");
		push @t, (0,$y,$y,"Sigma0");
		push @t, (1,1-$y,1-$y,"U0");
		push @t, (1/2,1/2+$z,$z,"A0");
		push @t, (1/2,1/2-$z,1-$z,"C0");
		push @t, (1/2,1/2,1/2,"L");
		if ($inversion ne "Y"){
			push @t, (-1,-1/2,-1/2,"T'");
			push @t, (-1/2,-1/2,0,"Z'");
			push @t, (-1/2,0,-1/2,"Y'");
			push @t, (0,-$y,-$y,"Sigma0'");
			push @t, (-1,-1+$y,-1+$y,"U0'");
			push @t, (-1/2,-1/2-$z,-$z,"A0'");
			push @t, (-1/2,-1/2+$z,-1+$z,"C0'");
			push @t, (-1/2,-1/2,-1/2,"L'");
		}
	} elsif ($latticetype eq "oF2"){
#oF2:c^-2>a^-2+b^-2
		my $z=(1+$c*$c/$a/$a-$c*$c/$b/$b)/4;
		my $y=(1+$c*$c/$a/$a+$c*$c/$b/$b)/4;
		push @t, (0,1/2,1/2,"T");
		push @t, (1/2,1/2,1,"Z");
		push @t, (1/2,0,1/2,"Y");
		push @t, ($y,$y,0,"Lambda0");
		push @t, (1-$y,1-$y,1,"Q0");
		push @t, (1/2-$z,1-$z,1/2,"G0");
		push @t, (1/2+$z,$z,1/2,"H0");
		push @t, (1/2,1/2,1/2,"L");
		if ($inversion ne "Y"){
			push @t, (0,-1/2,-1/2,"T'");
			push @t, (-1/2,-1/2,-1,"Z'");
			push @t, (-1/2,0,-1/2,"Y'");
			push @t, (-$y,-$y,0,"Lambda0'");
			push @t, (-1+$y,-1+$y,-1,"Q0'");
			push @t, (-1/2+$z,-1+$z,-1/2,"G0'");
			push @t, (-1/2-$z,-$z,-1/2,"H0'");
			push @t, (-1/2,-1/2,-1/2,"L'");
		}
	} elsif ($latticetype eq "oF3"){
#oF3:a^-2, b^-2, c^-2 edges of triangle
		my $y=(1+$a*$a/$b/$b-$a*$a/$c/$c)/4;
		my $d=(1+$b*$b/$a/$a-$b*$b/$c/$c)/4;
		my $p=(1+$c*$c/$b/$b-$c*$c/$a/$a)/4;
		push @t, (0,1/2,1/2,"T");
		push @t, (1/2,1/2,0,"Z");
		push @t, (1/2,0,1/2,"Y");
		push @t, (1/2,1/2+$y,$y,"A0");
		push @t, (1/2,1/2-$y,1-$y,"C0");
		push @t, (1/2+$d,1/2,$d,"B0");
		push @t, (1/2-$d,1/2,1-$d,"D0");
		push @t, ($p,1/2+$p,1/2,"G0");
		push @t, (1-$p,1/2-$p,1/2,"H0");
		push @t, (1/2,1/2,1/2,"L");
		if ($inversion ne "Y"){
			push @t, (0,-1/2,-1/2,"T'");
			push @t, (-1/2,-1/2,0,"Z'");
			push @t, (-1/2,0,-1/2,"Y'");
			push @t, (-1/2,-1/2-$y,-$y,"A0'");
			push @t, (-1/2,-1/2+$y,-1+$y,"C0'");
			push @t, (-1/2-$d,-1/2,-$d,"B0'");
			push @t, (-1/2+$d,-1/2,-1+$d,"D0'");
			push @t, (-$p,-1/2-$p,-1/2,"G0'");
			push @t, (-1+$p,-1/2+$p,-1/2,"H0'");
			push @t, (-1/2,-1/2,-1/2,"L'");
		}
	} elsif ($latticetype eq "oI1"){
#oI1:c largest
		my $z=(1+$a*$a/$c/$c)/4;
		my $y=(1+$b*$b/$c/$c)/4;
		my $d=($b*$b-$a*$a)/$c/$c/4;
		my $m=($b*$b+$a*$a)/$c/$c/4;
		push @t, (1/2,1/2,-1/2,"X");
		push @t, (1/2,0,0,"S");
		push @t, (0,1/2,0,"R");
		push @t, (0,0,1/2,"T");
		push @t, (1/4,1/4,1/4,"W");
		push @t, (-$z,$z,$z,"Sigma0");
		push @t, ($z,1-$z,-$z,"F2");
		push @t, ($y,-$y,$y,"Y0");
		push @t, (1-$y,$y,-$y,"U0");
		push @t, (-$m,$m,1/2-$d,"L0");
		push @t, ($m,-$m,1/2+$d,"M0");
		push @t, (1/2-$d,1/2+$d,-$m,"J0");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,1/2,"X'");
			push @t, (-1/2,0,0,"S'");
			push @t, (0,-1/2,0,"R'");
			push @t, (0,0,-1/2,"T'");
			push @t, (-1/4,-1/4,-1/4,"W'");
			push @t, ($z,-$z,-$z,"Sigma0'");
			push @t, (-$z,-1+$z,$z,"F2'");
			push @t, (-$y,$y,-$y,"Y0'");
			push @t, (-1+$y,-$y,$y,"U0'");
			push @t, ($m,-$m,-1/2+$d,"L0'");
			push @t, (-$m,$m,-1/2-$d,"M0'");
			push @t, (-1/2+$d,-1/2-$d,$m,"J0'");
		}
	} elsif ($latticetype eq "oI2"){
#oI2:a largest
		my $z=(1+$b*$b/$a/$a)/4;
		my $y=(1+$c*$c/$a/$a)/4;
		my $d=($c*$c-$b*$b)/$a/$a/4;
		my $m=($c*$c+$b*$b)/$a/$a/4;
		push @t, (-1/2,1/2,1/2,"X");
		push @t, (1/2,0,0,"S");
		push @t, (0,1/2,0,"R");
		push @t, (0,0,1/2,"T");
		push @t, (1/4,1/4,1/4,"W");
		push @t, ($z,-$z,$z,"Y0");
		push @t, (-$z,$z,1-$z,"U2");
		push @t, ($y,$y,-$y,"Lambda0");
		push @t, (-$y,1-$y,$y,"G2");
		push @t, (1/2-$d,-$m,$m,"K");
		push @t, (1/2+$d,$m,-$m,"K2");
		push @t, (-$m,1/2-$d,1/2+$d,"K4");
		if ($inversion ne "Y"){
			push @t, (1/2,-1/2,-1/2,"X'");
			push @t, (-1/2,0,0,"S'");
			push @t, (0,-1/2,0,"R'");
			push @t, (0,0,-1/2,"T'");
			push @t, (-1/4,-1/4,-1/4,"W'");
			push @t, (-$z,$z,-$z,"Y0'");
			push @t, ($z,-$z,-1+$z,"U2'");
			push @t, (-$y,-$y,$y,"Lambda0'");
			push @t, ($y,-1+$y,-$y,"G2'");
			push @t, (-1/2+$d,$m,-$m,"K'");
			push @t, (-1/2-$d,-$m,$m,"K2'");
			push @t, ($m,-1/2+$d,-1/2-$d,"K4'");
		}
	} elsif ($latticetype eq "oI3"){
#oI3:b largest
		my $z=(1+$c*$c/$b/$b)/4;
		my $y=(1+$a*$a/$b/$b)/4;
		my $d=($a*$a-$c*$c)/$b/$b/4;
		my $m=($a*$a+$c*$c)/$b/$b/4;
		push @t, (1/2,-1/2,1/2,"X");
		push @t, (1/2,0,0,"S");
		push @t, (0,1/2,0,"R");
		push @t, (0,0,1/2,"T");
		push @t, (1/4,1/4,1/4,"W");
		push @t, (-$y,$y,$y,"Sigma0");
		push @t, ($y,-$y,1-$y,"F0");
		push @t, ($z,$z,-$z,"Lambda0");
		push @t, (1-$z,-$z,$z,"G0");
		push @t, ($m,1/2-$d,-$m,"V0");
		push @t, (-$m,1/2+$d,$m,"H0");
		push @t, (1/2+$d,-$m,1/2-$d,"H2");
		if ($inversion ne "Y"){
			push @t, (-1/2,1/2,-1/2,"X'");
			push @t, (-1/2,0,0,"S'");
			push @t, (0,-1/2,0,"R'");
			push @t, (0,0,-1/2,"T'");
			push @t, (-1/4,-1/4,-1/4,"W'");
			push @t, ($y,-$y,-$y,"Sigma0'");
			push @t, (-$y,$y,-1+$y,"F0'");
			push @t, (-$z,-$z,$z,"Lambda0'");
			push @t, (-1+$z,$z,-$z,"G0'");
			push @t, (-$m,-1/2+$d,$m,"V0'");
			push @t, ($m,-1/2-$d,-$m,"H0'");
			push @t, (-1/2-$d,$m,-1/2+$d,"H2'");
		}
	} elsif (($latticetype eq "oC1") || ($latticetype eq "oA1")){
#oC1:a<b oA1:b<c
		my $z=(1+$a*$a/$b/$b)/4;
		$z=(1+$b*$b/$c/$c)/4 if ($latticetype eq "oA1");
		push @t, (-1/2,1/2,0,"Y");
		push @t, (-1/2,1/2,1/2,"T");
		push @t, (0,0,1/2,"Z");
		push @t, (0,1/2,0,"S");
		push @t, (0,1/2,1/2,"R");
		push @t, ($z,$z,0,"Sigma0");
		push @t, (-$z,1-$z,0,"C0");
		push @t, ($z,$z,1/2,"A0");
		push @t, (-$z,1-$z,1/2,"E0");
		if ($inversion ne "Y"){
			push @t, (1/2,-1/2,0,"Y'");
			push @t, (1/2,-1/2,-1/2,"T'");
			push @t, (0,0,-1/2,"Z'");
			push @t, (0,-1/2,0,"S'");
			push @t, (0,-1/2,-1/2,"R'");
			push @t, (-$z,-$z,0,"Sigma0'");
			push @t, ($z,-1+$z,0,"C0'");
			push @t, (-$z,-$z,-1/2,"A0'");
			push @t, ($z,-1+$z,-1/2,"E0'");
		}
	} elsif (($latticetype eq "oC2") || ($latticetype eq "oA2")){
#oC2:a>b oA2:b>c
		my $z=(1+$b*$b/$a/$a)/4;
		$z=(1+$c*$c/$b/$b)/4 if ($latticetype eq "oA2");
		push @t, (1/2,1/2,0,"Y");
		push @t, (1/2,1/2,1/2,"T");
		push @t, (1/2,1/2,-1/2,"T2");
		push @t, (0,0,1/2,"Z");
		push @t, (0,0,-1/2,"Z2");
		push @t, (0,1/2,0,"S");
		push @t, (0,1/2,1/2,"R");
		push @t, (0,1/2,-1/2,"R2");
		push @t, (-$z,$z,0,"Delta0");
		push @t, ($z,1-$z,0,"F0");
		push @t, (-$z,$z,1/2,"B0");
		push @t, (-$z,$z,-1/2,"B2");
		push @t, ($z,1-$z,1/2,"G0");
		push @t, ($z,1-$z,-1/2,"G2");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,0,"Y'");
			push @t, (-1/2,-1/2,-1/2,"T'");
			push @t, (-1/2,-1/2,1/2,"T2'");
			push @t, (0,0,-1/2,"Z'");
			push @t, (0,0,1/2,"Z2'");
			push @t, (0,-1/2,0,"S'");
			push @t, (0,-1/2,-1/2,"R'");
			push @t, (0,-1/2,1/2,"R2'");
			push @t, ($z,-$z,0,"Delta0'");
			push @t, (-$z,-1+$z,0,"F0'");
			push @t, ($z,-$z,-1/2,"B0'");
			push @t, ($z,-$z,1/2,"B2'");
			push @t, (-$z,-1+$z,-1/2,"G0'");
			push @t, (-$z,-1+$z,1/2,"G2'");
		}
	} elsif ($bravais eq "hP"){
		push @t, (0,0,1/2,"A");
		push @t, (1/3,1/3,0,"K");
		push @t, (1/3,1/3,1/2,"H");
		push @t, (1/3,1/3,-1/2,"H2");
		push @t, (1/2,0,0,"M");	
		push @t, (1/2,0,1/2,"L");
		if ($inversion ne "Y"){
			push @t, (0,0,-1/2,"A'");
			push @t, (-1/3,-1/3,0,"K'");
			push @t, (-1/3,-1/3,-1/2,"H'");
			push @t, (-1/3,-1/3,1/2,"H2'");
			push @t, (-1/2,0,0,"M'");	
			push @t, (-1/2,0,-1/2,"L'");
		}
	} elsif ($latticetype eq "hR1"){
#hR1:sqrt(3)a<sqrt(2)c
		my $d=$a*$a/$c/$c/4;
		my $y=5/6-2*$d;
		my $n=1/3+$d;
		push @t, (1/2,1/2,1/2,"T");
		push @t, (1/2,0,0,"L");
		push @t, (0,-1/2,0,"L2");
		push @t, (0,0,-1/2,"L4");
		push @t, (1/2,0,1/2,"F");
		push @t, (1/2,1/2,0,"F2");
		push @t, ($n,-$n,0,"S0");
		push @t, (1-$n,0,$n,"S2");
		push @t, ($n,0,-$n,"S4");
		push @t, (1-$n,$n,0,"S6");
		push @t, (1/2,-1+$y,1-$y,"H0");
		push @t, ($y,1-$y,1/2,"H2");
		push @t, ($y,1/2,1-$y,"H4");
		push @t, (1/2,1-$y,-1+$y,"H6");
		push @t, ($n,-1+$y,$n,"M");
		push @t, (1-$n,1-$y,1-$n,"M2");
		push @t, ($y,$n,$n,"M4");
		push @t, (1-$n,1-$n,1-$y,"M6");
		push @t, ($n,$n,-1+$y,"M8");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,-1/2,"T'");
			push @t, (-1/2,0,0,"L'");
			push @t, (0,1/2,0,"L2'");
			push @t, (0,0,1/2,"L4'");
			push @t, (-1/2,0,-1/2,"F'");
			push @t, (-1/2,-1/2,0,"F2'");
			push @t, (-$n,$n,0,"S0'");
			push @t, (-1+$n,0,-$n,"S2'");
			push @t, (-$n,0,$n,"S4'");
			push @t, (-1+$n,-$n,0,"S6'");
			push @t, (-1/2,1-$y,-1+$y,"H0'");
			push @t, (-$y,-1+$y,-1/2,"H2'");
			push @t, (-$y,-1/2,-1+$y,"H4'");
			push @t, (-1/2,-1+$y,1-$y,"H6'");
			push @t, (-$n,1-$y,-$n,"M'");
			push @t, (-1+$n,-1+$y,-1+$n,"M2'");
			push @t, (-$y,-$n,-$n,"M4'");
			push @t, (-1+$n,-1+$n,-1+$y,"M6'");
			push @t, (-$n,-$n,1-$y,"M8'");
		}
	} elsif ($latticetype eq "hR2"){
#hR2:sqrt(3)a>sqrt(2)c
		my $z=1/6-$c*$c/$a/$a/9;
		my $y=1/2-2*$z;
		my $n=1/2+$z;
		push @t, (1/2,-1/2,1/2,"T");
		push @t, ($y,-1+$y,$y,"P0");
		push @t, ($y,$y,$y,"P2");
		push @t, (1-$y,-$y,-$y,"R0");
		push @t, (1-$n,-$n,1-$n,"M");
		push @t, ($n,-1+$n,-1+$n,"M2");
		push @t, (1/2,0,0,"L");
		push @t, (1/2,-1/2,0,"F");
		if ($inversion ne "Y"){
			push @t, (-1/2,1/2,-1/2,"T'");
			push @t, (-$y,1-$y,-$y,"P0'");
			push @t, (-$y,-$y,-$y,"P2'");
			push @t, (-1+$y,$y,$y,"R0'");
			push @t, (-1+$n,$n,-1+$n,"M'");
			push @t, (-$n,1-$n,1-$n,"M2'");
			push @t, (-1/2,0,0,"L'");
			push @t, (-1/2,1/2,0,"F'");
		}
	} elsif ($bravais eq "mP"){
		my $y=(1+$a*$cbe/$c)/(2*$sbe*$sbe);
		my $n=1/2+$y*$c*$cbe/$a;
		push @t, (0,1/2,0,"Z");
		push @t, (0,0,1/2,"B");
		push @t, (0,0,-1/2,"B2");
		push @t, (1/2,0,0,"Y");
		push @t, (-1/2,0,0,"Y2");
		push @t, (1/2,1/2,0,"C");
		push @t, (-1/2,1/2,0,"C2");
		push @t, (0,1/2,1/2,"D");
		push @t, (0,1/2,-1/2,"D2");
		push @t, (-1/2,0,1/2,"A");
		push @t, (-1/2,1/2,1/2,"E");
		push @t, (-$y,0,1-$n,"H0");
		push @t, (-1+$y,0,$n,"H2");
		push @t, (-$y,0,-$n,"H4");
		push @t, (-$y,1/2,1-$n,"M");
		push @t, (-1+$y,1/2,$n,"M2");
		push @t, (-$y,1/2,-$n,"M4");
		if ($inversion ne "Y"){
			push @t, (0,-1/2,0,"Z'");
			push @t, (0,0,-1/2,"B'");
			push @t, (0,0,1/2,"B2'");
			push @t, (-1/2,0,0,"Y'");
			push @t, (1/2,0,0,"Y2'");
			push @t, (-1/2,-1/2,0,"C'");
			push @t, (1/2,-1/2,0,"C2'");
			push @t, (0,-1/2,-1/2,"D'");
			push @t, (0,-1/2,1/2,"D2'");
			push @t, (1/2,0,-1/2,"A'");
			push @t, (1/2,-1/2,-1/2,"E'");
			push @t, ($y,0,-1+$n,"H'");
			push @t, (1-$y,0,-$n,"H2'");
			push @t, ($y,0,$n,"H4'");
			push @t, ($y,-1/2,-1+$n,"M'");
			push @t, (1-$y,-1/2,-$n,"M2'");
			push @t, ($y,-1/2,$n,"M4'");
		}
	} elsif ($latticetype eq "mC1"){
#mC1:b<asin(beta)
		my $z=(2+$a*$cbe/$c)/4/$sbe/$sbe;
		my $y=1/2-2*$z*$c*$cbe/$a;
		my $s=3/4-$b*$b/4/$a/$a/$sbe/$sbe;
		my $p=$s-(3/4-$s)*$a*$cbe/$c;
		push @t, (-1/2,1/2,0,"Y2");
		push @t, (1/2,-1/2,0,"Y4");
		push @t, (0,0,1/2,"A");
		push @t, (-1/2,1/2,1/2,"M2");
		push @t, (1/2,0,0,"V");
		push @t, (0,1/2,0,"V2");
		push @t, (0,1/2,1/2,"L2");
		push @t, (1-$s,1-$s,0,"C");
		push @t, (-1+$s,$s,0,"C2");
		push @t, ($s,-1+$s,0,"C4");
		push @t, (-1+$p,$p,1/2,"D");
		push @t, (1-$p,1-$p,1/2,"D2");
		push @t, (-1+$z,1-$z,1-$y,"E");
		push @t, (-$z,$z,$y,"E2");
		push @t, ($z,-$z,1-$y,"E4");
		if ($inversion ne "Y"){
			push @t, (1/2,-1/2,0,"Y2'");
			push @t, (-1/2,1/2,0,"Y4'");
			push @t, (0,0,-1/2,"A'");
			push @t, (1/2,-1/2,-1/2,"M2'");
			push @t, (-1/2,0,0,"V'");
			push @t, (0,-1/2,0,"V2'");
			push @t, (0,-1/2,-1/2,"L2'");
			push @t, (-1+$s,-1+$s,0,"C'");
			push @t, (1-$s,-$s,0,"C2'");
			push @t, (-$s,1-$s,0,"C4'");
			push @t, (1-$p,-$p,-1/2,"D'");
			push @t, (-1+$p,-1+$p,-1/2,"D2'");
			push @t, (1-$z,-1+$z,-1+$y,"E'");
			push @t, ($z,-$z,-$y,"E2'");
			push @t, (-$z,$z,-1+$y,"E4'");
		}
	} elsif ($latticetype eq "mC2"){
#mC2:b>asin(beta) BZ 12-face
		my $m=(1+$a*$a/$b/$b)/4;
		my $d=-$a*$c*$cbe/2/$b/$b;
		my $z=$a*$a/$b/$b/4+(1+$a*$cbe/$c)/4/$sbe/$sbe;
		my $y=1/2-2*$z*$c*$cbe/$a;
		my $p=1+$z-2*$m;
		my $s=$y-2*$d;
		push @t, (1/2,1/2,0,"Y");
		push @t, (0,0,1/2,"A");
		push @t, (1/2,1/2,1/2,"M");
		push @t, (0,1/2,0,"V2");
		push @t, (0,1/2,1/2,"L2");
		push @t, (-1+$p,1-$p,1-$s,"F");
		push @t, (1-$p,$p,$s,"F2");
		push @t, ($p,1-$p,1-$s,"F4");
		push @t, (-$z,$z,$y,"H");
		push @t, ($z,1-$z,1-$y,"H2");
		push @t, ($z,-$z,1-$y,"H4");
		push @t, (-$m,$m,$d,"G");
		push @t, ($m,1-$m,-$d,"G2");
		push @t, ($m,-$m,-$d,"G4");
		push @t, (1-$m,$m,$d,"G6");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,0,"Y'");
			push @t, (0,0,-1/2,"A'");
			push @t, (-1/2,-1/2,-1/2,"M'");
			push @t, (0,-1/2,0,"V2'");
			push @t, (0,-1/2,-1/2,"L2'");
			push @t, (1-$p,-1+$p,-1+$s,"F'");
			push @t, (-1+$p,-$p,-$s,"F2'");
			push @t, (-$p,-1+$p,-1+$s,"F4'");
			push @t, ($z,-$z,-$y,"H'");
			push @t, (-$z,-1+$z,-1+$y,"H2'");
			push @t, (-$z,$z,-1+$y,"H4'");
			push @t, ($m,-$m,-$d,"G'");
			push @t, (-$m,-1+$m,$d,"G2'");
			push @t, (-$m,$m,$d,"G4'");
			push @t, (-1+$m,-$m,-$d,"G6'");
		}
	} elsif ($latticetype eq "mC3"){
#mC3:b>asin(beta) BZ 14-face
		my $z=($a*$a/$b/$b+(1+$a*$cbe/$c)/$sbe/$sbe)/4;
		my $y=1/2-2*$z*$c*$cbe/$a;
		my $m=$y/2+$a*$a/4/$b/$b+$a*$c*$cbe/2/$b/$b;
		my $n=2*$m-$z;
		my $o=(1-4*$n+$a*$a*$sbe*$sbe/$b/$b)*$c/(2*$a*$cbe);
		my $d=$o/2-1/4-$z*$c*$cbe/$a;
		my $r=1-$z*$b*$b/$a/$a;
		push @t, (1/2,1/2,0,"Y");
		push @t, (0,0,1/2,"A");
		push @t, (-1/2,1/2,1/2,"M2");
		push @t, (1/2,0,0,"V");
		push @t, (0,1/2,0,"V2");
		push @t, (0,1/2,1/2,"L2");
		push @t, (-1+$r,$r,1/2,"I");
		push @t, (1-$r,1-$r,1/2,"I2");
		push @t, (-$n,$n,$o,"K");
		push @t, (-1+$n,1-$n,1-$o,"K2");
		push @t, (1-$n,$n,$o,"K4");
		push @t, (-$z,$z,$y,"H");
		push @t, ($z,1-$z,1-$y,"H2");
		push @t, ($z,-$z,1-$y,"H4");
		push @t, (-$m,$m,$d,"N");
		push @t, ($m,1-$m,-$d,"N2");
		push @t, ($m,-$m,-$d,"N4");
		push @t, (1-$m,$m,$d,"N6");
		if ($inversion ne "Y"){
			push @t, (-1/2,-1/2,0,"Y'");
			push @t, (0,0,-1/2,"A'");
			push @t, (1/2,-1/2,-1/2,"M2'");
			push @t, (-1/2,0,0,"V'");
			push @t, (0,-1/2,0,"V2'");
			push @t, (0,-1/2,-1/2,"L2'");
			push @t, (1-$r,-$r,-1/2,"I'");
			push @t, (-1+$r,-1+$r,-1/2,"I2'");
			push @t, ($n,-$n,-$o,"K'");
			push @t, (1-$n,-1+$n,-1+$o,"K2'");
			push @t, (-1+$n,-$n,-$o,"K4'");
			push @t, ($z,-$z,-$y,"H'");
			push @t, (-$z,-1+$z,-1+$y,"H2'");
			push @t, (-$z,$z,-1+$y,"H4'");
			push @t, ($m,-$m,-$d,"N'");
			push @t, (-$m,-1+$m,$d,"N2'");
			push @t, (-$m,$m,$d,"N4'");
			push @t, (-1+$m,-$m,-$d,"N6'");
		}
	} elsif ($bravais eq "aP"){
#triclinic: REDUCED NI SURU!
		my $type=&get_hinuma_reduced;
		if ($type eq "Obtuse"){
			push @t, (0,0,1/2,"Z");
			push @t, (0,1/2,0,"Y");
			push @t, (1/2,0,0,"X");
			push @t, (1/2,1/2,0,"V");
			push @t, (1/2,0,1/2,"U");
			push @t, (1/2,1/2,1/2,"T");
			push @t, (0,1/2,1/2,"R");
			if ($inversion ne "Y"){
				push @t, (0,0,-1/2,"Z'");
				push @t, (0,-1/2,0,"Y'");
				push @t, (-1/2,0,0,"X'");
				push @t, (-1/2,-1/2,0,"V'");
				push @t, (-1/2,0,-1/2,"U'");
				push @t, (-1/2,-1/2,-1/2,"T'");
				push @t, (0,-1/2,-1/2,"R'");
			}
		} elsif ($type eq "Acute"){
			push @t, (0,0,1/2,"Z");
			push @t, (0,1/2,0,"Y");
			push @t, (0,-1/2,0,"Y2");
			push @t, (1/2,0,0,"X");
			push @t, (1/2,-1/2,0,"V2");
			push @t, (-1/2,0,1/2,"U2");
			push @t, (0,-1/2,1/2,"T2");
			push @t, (-1/2,-1/2,1/2,"R2");
			if ($inversion ne "Y"){
				push @t, (0,0,-1/2,"Z'");
				push @t, (0,-1/2,0,"Y'");
				push @t, (0,1/2,0,"Y2'");
				push @t, (-1/2,0,0,"X'");
				push @t, (-1/2,1/2,0,"V2'");
				push @t, (1/2,0,-1/2,"U2'");
				push @t, (0,1/2,-1/2,"T2'");
				push @t, (1/2,1/2,-1/2,"R2'");
			}
		} else {
			die ("bz_point: triclinic no type ga okasii \n");
		}
	} else {
		die ("bz_point: lattice type ga okasii \n");
	}
	&get_crystallographic_cell("primitive") if ($centering ne "P");
	return (@t);
}

#Kessyougaku no kpath wo motomeru
#kpf huu ni tsubo_temp_get_kpath ni syuturyoku
#file mei ha betu no basyo de kaeru!
sub get_kpath_kpf{
#Conventional wo nyuusyu: k-point data wo motomeru
	my @points=&bz_point;
	my $extended_bravais=shift(@points);
	my $inversion=substr($extended_bravais,3,3);
	my $bravais=substr($extended_bravais,0,2);
#	my $centering=substr($extended_bravais,1,1);
	my $latticetype=substr($extended_bravais,0,3);
	open OUT, ">tsubo_temp_get_kpath";
	print OUT ("Real form of k-point coordinates (kx,ky,kz,label):\n");
	$inversion="Y" if ($_[0] eq "Y");
	&get_kpath_kpf_print(1,0,"-",@points);
	if ($latticetype eq "cP1"){
#cP1: 0-3-2-0-1-3|1-2(-4)
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,3,"|",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		&get_kpath_kpf_print(8,2,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(9,4,"X",@points);
		} else {
			&get_kpath_kpf_print(9,4,"|",@points);
			&get_kpath_kpf_print(10,0,"-",@points);
			&get_kpath_kpf_print(11,7,"-",@points);
			&get_kpath_kpf_print(12,6,"-",@points);
			&get_kpath_kpf_print(13,0,"-",@points);
			&get_kpath_kpf_print(14,5,"-",@points);
			&get_kpath_kpf_print(15,7,"|",@points);
			&get_kpath_kpf_print(16,5,"-",@points);
			&get_kpath_kpf_print(17,6,"-",@points);
			&get_kpath_kpf_print(18,8,"X",@points);
		}
	} elsif ($latticetype eq "cP2"){
#cP2: 0-3-2-0-1-3|1-2
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,3,"|",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(8,2,"X",@points);
		} else {
			&get_kpath_kpf_print(8,2,"|",@points);
			&get_kpath_kpf_print(9,0,"-",@points);
			&get_kpath_kpf_print(10,7,"-",@points);
			&get_kpath_kpf_print(11,6,"-",@points);
			&get_kpath_kpf_print(12,0,"-",@points);
			&get_kpath_kpf_print(13,5,"-",@points);
			&get_kpath_kpf_print(14,7,"|",@points);
			&get_kpath_kpf_print(15,5,"-",@points);
			&get_kpath_kpf_print(16,6,"X",@points);

		}
	} elsif ($latticetype eq "cF1"){
#cF1:0-1-6|5-0-2-3-1-4
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,6,"|",@points);
		&get_kpath_kpf_print(4,5,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,2,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,1,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(9,4,"X",@points);
		} else {
			&get_kpath_kpf_print(9,4,"|",@points);
			&get_kpath_kpf_print(10,0,"-",@points);
			&get_kpath_kpf_print(11,7,"-",@points);
			&get_kpath_kpf_print(12,12,"|",@points);
			&get_kpath_kpf_print(13,11,"-",@points);
			&get_kpath_kpf_print(14,0,"-",@points);
			&get_kpath_kpf_print(15,8,"-",@points);
			&get_kpath_kpf_print(16,9,"-",@points);
			&get_kpath_kpf_print(17,7,"-",@points);
			&get_kpath_kpf_print(18,10,"X",@points);
		}
	} elsif ($latticetype eq "cF2"){
#cF2:0-1-6|5-0-2-3-1
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,6,"|",@points);
		&get_kpath_kpf_print(4,5,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,2,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(8,1,"X",@points);
		} else {
			&get_kpath_kpf_print(8,1,"|",@points);
			&get_kpath_kpf_print(9,0,"-",@points);
			&get_kpath_kpf_print(10,7,"-",@points);
			&get_kpath_kpf_print(11,12,"|",@points);
			&get_kpath_kpf_print(12,11,"-",@points);
			&get_kpath_kpf_print(13,0,"-",@points);
			&get_kpath_kpf_print(14,8,"-",@points);
			&get_kpath_kpf_print(15,9,"-",@points);
			&get_kpath_kpf_print(16,7,"X",@points);
		}
	} elsif ($bravais eq "cI"){
#cI: 0-1-3-0-2-1|2-3
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,2,"-",@points);
		&get_kpath_kpf_print(6,1,"|",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(8,3,"X",@points);
		} else {
			&get_kpath_kpf_print(8,3,"X",@points);
			&get_kpath_kpf_print(9,0,"-",@points);
			&get_kpath_kpf_print(10,4,"-",@points);
			&get_kpath_kpf_print(11,6,"-",@points);
			&get_kpath_kpf_print(12,0,"-",@points);
			&get_kpath_kpf_print(13,5,"-",@points);
			&get_kpath_kpf_print(14,4,"|",@points);
			&get_kpath_kpf_print(15,5,"-",@points);
			&get_kpath_kpf_print(15,6,"X",@points);
		}
	} elsif ($bravais eq "tP"){
#tP:0-5-2-0-1-4-3-1|5-4|2-3
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,4,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,5,"-",@points);
		&get_kpath_kpf_print(10,4,"|",@points);
		&get_kpath_kpf_print(11,2,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(12,3,"X",@points);
		} else {
			&get_kpath_kpf_print(12,3,"|",@points);
			&get_kpath_kpf_print(13,0,"-",@points);
			&get_kpath_kpf_print(14,10,"-",@points);
			&get_kpath_kpf_print(15,7,"-",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,6,"-",@points);
			&get_kpath_kpf_print(18,9,"-",@points);
			&get_kpath_kpf_print(19,8,"-",@points);
			&get_kpath_kpf_print(20,6,"|",@points);
			&get_kpath_kpf_print(21,10,"-",@points);
			&get_kpath_kpf_print(22,9,"-",@points);
			&get_kpath_kpf_print(23,7,"-",@points);
			&get_kpath_kpf_print(24,8,"X",@points);
		}
	} elsif ($latticetype eq "tI1"){
# tI1:c<a 0-2-1-0-4|5-1|2-3-6-0
		&get_kpath_kpf_print(2,2,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,4,"|",@points);
		&get_kpath_kpf_print(6,5,"-",@points);
		&get_kpath_kpf_print(7,1,"|",@points);
		&get_kpath_kpf_print(8,2,"-",@points);
		&get_kpath_kpf_print(9,3,"-",@points);
		&get_kpath_kpf_print(10,6,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(11,0,"X",@points);
		} else {
			&get_kpath_kpf_print(11,0,"-",@points);
			&get_kpath_kpf_print(12,8,"-",@points);
			&get_kpath_kpf_print(13,7,"-",@points);
			&get_kpath_kpf_print(14,0,"-",@points);
			&get_kpath_kpf_print(15,10,"|",@points);
			&get_kpath_kpf_print(16,11,"-",@points);
			&get_kpath_kpf_print(17,7,"|",@points);
			&get_kpath_kpf_print(18,8,"-",@points);
			&get_kpath_kpf_print(19,9,"-",@points);
			&get_kpath_kpf_print(20,12,"-",@points);
			&get_kpath_kpf_print(21,0,"X",@points);
		}
	} elsif ($latticetype eq "tI2"){
# tI2:c>a 0-2-3-4-0-1-6|5-0|2-7|8-1
		&get_kpath_kpf_print(2,2,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,1,"-",@points);
		&get_kpath_kpf_print(7,6,"|",@points);
		&get_kpath_kpf_print(8,5,"-",@points);
		&get_kpath_kpf_print(9,0,"|",@points);
		&get_kpath_kpf_print(10,2,"-",@points);
		&get_kpath_kpf_print(11,7,"|",@points);
		&get_kpath_kpf_print(12,8,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(13,1,"X",@points);
		} else {
			&get_kpath_kpf_print(13,1,"|",@points);
			&get_kpath_kpf_print(14,0,"-",@points);
			&get_kpath_kpf_print(15,10,"-",@points);
			&get_kpath_kpf_print(16,11,"-",@points);
			&get_kpath_kpf_print(17,12,"-",@points);
			&get_kpath_kpf_print(18,0,"-",@points);
			&get_kpath_kpf_print(19,9,"-",@points);
			&get_kpath_kpf_print(20,14,"|",@points);
			&get_kpath_kpf_print(21,13,"-",@points);
			&get_kpath_kpf_print(22,0,"|",@points);
			&get_kpath_kpf_print(23,10,"-",@points);
			&get_kpath_kpf_print(24,15,"|",@points);
			&get_kpath_kpf_print(25,16,"-",@points);
			&get_kpath_kpf_print(26,9,"X",@points);
		}
	} elsif ($bravais eq "oP"){
#oP:0-1-5-4-0-2-3-7-6-2|1-3|4-6|5-7
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,5,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,2,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,7,"-",@points);
		&get_kpath_kpf_print(9,6,"-",@points);
		&get_kpath_kpf_print(10,2,"|",@points);
		&get_kpath_kpf_print(11,1,"-",@points);
		&get_kpath_kpf_print(12,3,"|",@points);
		&get_kpath_kpf_print(13,4,"-",@points);
		&get_kpath_kpf_print(14,6,"|",@points);
		&get_kpath_kpf_print(15,5,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(16,7,"X",@points);
		} else {
			&get_kpath_kpf_print(16,7,"|",@points);
			&get_kpath_kpf_print(17,0,"-",@points);
			&get_kpath_kpf_print(18,8,"-",@points);
			&get_kpath_kpf_print(19,12,"-",@points);
			&get_kpath_kpf_print(20,11,"-",@points);
			&get_kpath_kpf_print(21,0,"-",@points);
			&get_kpath_kpf_print(22,9,"-",@points);
			&get_kpath_kpf_print(23,10,"-",@points);
			&get_kpath_kpf_print(24,14,"-",@points);
			&get_kpath_kpf_print(25,13,"-",@points);
			&get_kpath_kpf_print(26,9,"|",@points);
			&get_kpath_kpf_print(27,8,"-",@points);
			&get_kpath_kpf_print(28,10,"|",@points);
			&get_kpath_kpf_print(29,11,"-",@points);
			&get_kpath_kpf_print(30,13,"|",@points);
			&get_kpath_kpf_print(31,12,"-",@points);
			&get_kpath_kpf_print(32,14,"X",@points);
		}
	} elsif ($latticetype eq "oF1"){
#oF1:a^-2>b^-2+c^-2  0-3-1-2-0-4|5-1|3-7|6-2|0-8
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,2,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,4,"|",@points);
		&get_kpath_kpf_print(7,5,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,3,"-",@points);
		&get_kpath_kpf_print(10,7,"|",@points);
		&get_kpath_kpf_print(11,6,"-",@points);
		&get_kpath_kpf_print(12,2,"|",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(14,8,"X",@points);
		} else {
			&get_kpath_kpf_print(14,8,"|",@points);
			&get_kpath_kpf_print(15,0,"-",@points);
			&get_kpath_kpf_print(16,11,"-",@points);
			&get_kpath_kpf_print(17,9,"-",@points);
			&get_kpath_kpf_print(18,10,"-",@points);
			&get_kpath_kpf_print(19,0,"-",@points);
			&get_kpath_kpf_print(20,12,"|",@points);
			&get_kpath_kpf_print(21,13,"-",@points);
			&get_kpath_kpf_print(22,9,"|",@points);
			&get_kpath_kpf_print(23,11,"-",@points);
			&get_kpath_kpf_print(24,15,"|",@points);
			&get_kpath_kpf_print(25,14,"-",@points);
			&get_kpath_kpf_print(26,10,"|",@points);
			&get_kpath_kpf_print(27,0,"-",@points);
			&get_kpath_kpf_print(28,16,"X",@points);
		}
	} elsif ($latticetype eq "oF2"){
#oF2:c^-2>a^-2+b^-2 0-1-2-3-0-4|5-2|1-6|7-3|0-8
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,3,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,4,"|",@points);
		&get_kpath_kpf_print(7,5,"-",@points);
		&get_kpath_kpf_print(8,2,"|",@points);
		&get_kpath_kpf_print(9,1,"-",@points);
		&get_kpath_kpf_print(10,6,"|",@points);
		&get_kpath_kpf_print(11,7,"-",@points);
		&get_kpath_kpf_print(12,3,"|",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(14,8,"X",@points);
		} else {
			&get_kpath_kpf_print(14,8,"|",@points);
			&get_kpath_kpf_print(15,0,"-",@points);
			&get_kpath_kpf_print(16,9,"-",@points);
			&get_kpath_kpf_print(17,10,"-",@points);
			&get_kpath_kpf_print(18,11,"-",@points);
			&get_kpath_kpf_print(19,0,"-",@points);
			&get_kpath_kpf_print(20,12,"|",@points);
			&get_kpath_kpf_print(21,13,"-",@points);
			&get_kpath_kpf_print(22,10,"|",@points);
			&get_kpath_kpf_print(23,9,"-",@points);
			&get_kpath_kpf_print(24,14,"|",@points);
			&get_kpath_kpf_print(25,15,"-",@points);
			&get_kpath_kpf_print(26,11,"|",@points);
			&get_kpath_kpf_print(27,0,"-",@points);
			&get_kpath_kpf_print(28,16,"X",@points);
		}
	} elsif ($latticetype eq "oF3"){
#oF3:a^-2, b^-2, c^-2 edges of triangle
# 0-3-5|4-2-6|7-1-8|9-3|1-0-2|0-10
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,5,"|",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,2,"-",@points);
		&get_kpath_kpf_print(6,6,"|",@points);
		&get_kpath_kpf_print(7,7,"-",@points);
		&get_kpath_kpf_print(8,1,"-",@points);
		&get_kpath_kpf_print(9,8,"|",@points);
		&get_kpath_kpf_print(10,9,"-",@points);
		&get_kpath_kpf_print(11,3,"|",@points);
		&get_kpath_kpf_print(12,1,"-",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		&get_kpath_kpf_print(14,2,"|",@points);
		&get_kpath_kpf_print(15,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(16,10,"X",@points);
		} else {
			&get_kpath_kpf_print(16,10,"|",@points);
			&get_kpath_kpf_print(17,0,"-",@points);
			&get_kpath_kpf_print(18,13,"-",@points);
			&get_kpath_kpf_print(19,15,"|",@points);
			&get_kpath_kpf_print(20,14,"-",@points);
			&get_kpath_kpf_print(21,12,"-",@points);
			&get_kpath_kpf_print(22,16,"|",@points);
			&get_kpath_kpf_print(23,17,"-",@points);
			&get_kpath_kpf_print(24,11,"-",@points);
			&get_kpath_kpf_print(25,18,"|",@points);
			&get_kpath_kpf_print(26,19,"-",@points);
			&get_kpath_kpf_print(27,13,"|",@points);
			&get_kpath_kpf_print(28,11,"-",@points);
			&get_kpath_kpf_print(29,0,"-",@points);
			&get_kpath_kpf_print(30,12,"|",@points);
			&get_kpath_kpf_print(31,0,"-",@points);
			&get_kpath_kpf_print(32,20,"X",@points);
		}
	} elsif ($bravais eq "oI"){
#oI:common 0-1-7|6-0-8|9-1|0-3-5-2-0-4-5
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,7,"|",@points);
		&get_kpath_kpf_print(4,6,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,8,"|",@points);
		&get_kpath_kpf_print(7,9,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,0,"-",@points);
		&get_kpath_kpf_print(10,3,"-",@points);
		&get_kpath_kpf_print(11,5,"-",@points);
		&get_kpath_kpf_print(12,2,"-",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		&get_kpath_kpf_print(14,4,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(15,5,"X",@points);
		} else {
			&get_kpath_kpf_print(15,5,"|",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,13,"-",@points);
			&get_kpath_kpf_print(18,19,"|",@points);
			&get_kpath_kpf_print(19,18,"-",@points);
			&get_kpath_kpf_print(20,0,"-",@points);
			&get_kpath_kpf_print(21,20,"|",@points);
			&get_kpath_kpf_print(22,21,"-",@points);
			&get_kpath_kpf_print(23,13,"|",@points);
			&get_kpath_kpf_print(24,0,"-",@points);
			&get_kpath_kpf_print(25,15,"-",@points);
			&get_kpath_kpf_print(26,17,"-",@points);
			&get_kpath_kpf_print(27,14,"-",@points);
			&get_kpath_kpf_print(28,0,"-",@points);
			&get_kpath_kpf_print(29,16,"-",@points);
			&get_kpath_kpf_print(30,17,"X",@points);
		}
	} elsif (($latticetype eq "oC1") || ($latticetype eq "oA1")){
#oC1:a<b oA1:b<c 0-1-7|6-0-3-8|9-2-1|0-4-5-3-2
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,7,"|",@points);
		&get_kpath_kpf_print(4,6,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,3,"-",@points);
		&get_kpath_kpf_print(7,8,"|",@points);
		&get_kpath_kpf_print(8,9,"-",@points);
		&get_kpath_kpf_print(9,2,"-",@points);
		&get_kpath_kpf_print(10,1,"|",@points);
		&get_kpath_kpf_print(11,0,"-",@points);
		&get_kpath_kpf_print(12,4,"-",@points);
		&get_kpath_kpf_print(13,5,"-",@points);
		&get_kpath_kpf_print(14,3,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(15,2,"X",@points);
		} else {
			&get_kpath_kpf_print(15,2,"|",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,10,"-",@points);
			&get_kpath_kpf_print(18,16,"|",@points);
			&get_kpath_kpf_print(19,15,"-",@points);
			&get_kpath_kpf_print(20,0,"-",@points);
			&get_kpath_kpf_print(21,12,"-",@points);
			&get_kpath_kpf_print(22,17,"|",@points);
			&get_kpath_kpf_print(23,18,"-",@points);
			&get_kpath_kpf_print(24,11,"-",@points);
			&get_kpath_kpf_print(25,10,"|",@points);
			&get_kpath_kpf_print(26,0,"-",@points);
			&get_kpath_kpf_print(27,13,"-",@points);
			&get_kpath_kpf_print(28,14,"-",@points);
			&get_kpath_kpf_print(29,12,"-",@points);
			&get_kpath_kpf_print(30,11,"-",@points);
		}
	} elsif (($latticetype eq "oC2") || ($latticetype eq "oA2")){
#oC2:a>b oA2:b>c 0-1-10|9-0-4-11|13-2-1|0-6-7-4-2
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,10,"|",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,4,"-",@points);
		&get_kpath_kpf_print(7,11,"|",@points);
		&get_kpath_kpf_print(8,13,"-",@points);
		&get_kpath_kpf_print(9,2,"-",@points);
		&get_kpath_kpf_print(10,1,"|",@points);
		&get_kpath_kpf_print(11,0,"-",@points);
		&get_kpath_kpf_print(12,6,"-",@points);
		&get_kpath_kpf_print(13,7,"-",@points);
		&get_kpath_kpf_print(14,4,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(15,2,"X",@points);
		} else {
			&get_kpath_kpf_print(15,2,"|",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,15,"-",@points);
			&get_kpath_kpf_print(18,24,"|",@points);
			&get_kpath_kpf_print(19,23,"-",@points);
			&get_kpath_kpf_print(20,0,"-",@points);
			&get_kpath_kpf_print(21,18,"-",@points);
			&get_kpath_kpf_print(22,25,"|",@points);
			&get_kpath_kpf_print(23,27,"-",@points);
			&get_kpath_kpf_print(24,16,"-",@points);
			&get_kpath_kpf_print(25,15,"|",@points);
			&get_kpath_kpf_print(26,0,"-",@points);
			&get_kpath_kpf_print(27,20,"-",@points);
			&get_kpath_kpf_print(28,21,"-",@points);
			&get_kpath_kpf_print(29,18,"-",@points);
			&get_kpath_kpf_print(30,16,"X",@points);
		}
	} elsif ($latticetype eq "hP1"){
#hP1 0-5-2-0-1-6-3-1|6-5|3-2-4
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,6,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,6,"-",@points);
		&get_kpath_kpf_print(10,5,"|",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		&get_kpath_kpf_print(12,2,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(13,4,"X",@points);
		} else {
			&get_kpath_kpf_print(13,4,"|",@points);
			&get_kpath_kpf_print(14,0,"-",@points);
			&get_kpath_kpf_print(15,11,"-",@points);
			&get_kpath_kpf_print(16,8,"-",@points);
			&get_kpath_kpf_print(17,0,"-",@points);
			&get_kpath_kpf_print(18,7,"-",@points);
			&get_kpath_kpf_print(19,12,"-",@points);
			&get_kpath_kpf_print(20,9,"-",@points);
			&get_kpath_kpf_print(21,7,"|",@points);
			&get_kpath_kpf_print(22,12,"-",@points);
			&get_kpath_kpf_print(23,11,"|",@points);
			&get_kpath_kpf_print(24,9,"-",@points);
			&get_kpath_kpf_print(25,8,"-",@points);
			&get_kpath_kpf_print(26,10,"X",@points);
		}
	} elsif ($latticetype eq "hP2"){
#hP2 0-5-2-0-1-6-3-1|6-5|3-2
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,6,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,6,"-",@points);
		&get_kpath_kpf_print(10,5,"|",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(12,2,"X",@points);
		} else {
			&get_kpath_kpf_print(12,2,"|",@points);
			&get_kpath_kpf_print(13,0,"-",@points);
			&get_kpath_kpf_print(14,11,"-",@points);
			&get_kpath_kpf_print(15,8,"-",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,7,"-",@points);
			&get_kpath_kpf_print(18,12,"-",@points);
			&get_kpath_kpf_print(19,9,"-",@points);
			&get_kpath_kpf_print(20,7,"|",@points);
			&get_kpath_kpf_print(21,12,"-",@points);
			&get_kpath_kpf_print(22,11,"|",@points);
			&get_kpath_kpf_print(23,9,"-",@points);
			&get_kpath_kpf_print(24,8,"X",@points);
		}
	} elsif ($latticetype eq "hR1"){
#hR1:sqrt(3)a<sqrt(2)c 0-1-12|11-2-0-7|8-5-0
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,12,"|",@points);
		&get_kpath_kpf_print(4,11,"-",@points);
		&get_kpath_kpf_print(5,2,"-",@points);
		&get_kpath_kpf_print(6,0,"-",@points);
		&get_kpath_kpf_print(7,7,"|",@points);
		&get_kpath_kpf_print(8,8,"-",@points);
		&get_kpath_kpf_print(9,5,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(10,0,"X",@points);
		} else {
			&get_kpath_kpf_print(10,0,"-",@points);
			&get_kpath_kpf_print(11,20,"-",@points);
			&get_kpath_kpf_print(12,31,"|",@points);
			&get_kpath_kpf_print(13,30,"-",@points);
			&get_kpath_kpf_print(14,21,"-",@points);
			&get_kpath_kpf_print(15,0,"-",@points);
			&get_kpath_kpf_print(16,26,"|",@points);
			&get_kpath_kpf_print(17,27,"-",@points);
			&get_kpath_kpf_print(18,24,"-",@points);
			&get_kpath_kpf_print(19,0,"X",@points);
		}
	} elsif ($latticetype eq "hR2"){
#hR2:sqrt(3)a>sqrt(2)c 0-7-1-2|3-0-8
		&get_kpath_kpf_print(2,7,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,2,"|",@points);
		&get_kpath_kpf_print(5,3,"-",@points);
		&get_kpath_kpf_print(6,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(7,8,"X",@points);
		} else {
			&get_kpath_kpf_print(7,8,"|",@points);
			&get_kpath_kpf_print(8,0,"-",@points);
			&get_kpath_kpf_print(9,15,"-",@points);
			&get_kpath_kpf_print(10,9,"-",@points);
			&get_kpath_kpf_print(11,10,"|",@points);
			&get_kpath_kpf_print(12,11,"-",@points);
			&get_kpath_kpf_print(13,0,"-",@points);
			&get_kpath_kpf_print(14,16,"X",@points);
		}
	} elsif ($bravais eq "mP"){
#mP: 0-1-8-2-0-10-11-1-7-5-0
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,8,"-",@points);
		&get_kpath_kpf_print(4,2,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,10,"-",@points);
		&get_kpath_kpf_print(7,11,"-",@points);
		&get_kpath_kpf_print(8,1,"-",@points);
		&get_kpath_kpf_print(9,7,"-",@points);
		&get_kpath_kpf_print(10,5,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(11,0,"X",@points);
		} else {
			&get_kpath_kpf_print(11,0,"-",@points);
			&get_kpath_kpf_print(12,18,"-",@points);
			&get_kpath_kpf_print(13,25,"-",@points);
			&get_kpath_kpf_print(14,19,"-",@points);
			&get_kpath_kpf_print(15,0,"-",@points);
			&get_kpath_kpf_print(16,27,"-",@points);
			&get_kpath_kpf_print(17,28,"-",@points);
			&get_kpath_kpf_print(18,18,"-",@points);
			&get_kpath_kpf_print(19,24,"-",@points);
			&get_kpath_kpf_print(20,22,"-",@points);
			&get_kpath_kpf_print(21,0,"X",@points);
		}
	} elsif ($latticetype eq "mC1"){
#mC1:b<asin(beta) 0-8|9-1-0-4-11|12-3-0|7-0-6
		&get_kpath_kpf_print(2,8,"|",@points);
		&get_kpath_kpf_print(3,9,"-",@points);
		&get_kpath_kpf_print(4,1,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,4,"-",@points);
		&get_kpath_kpf_print(7,11,"|",@points);
		&get_kpath_kpf_print(8,12,"-",@points);
		&get_kpath_kpf_print(9,3,"-",@points);
		&get_kpath_kpf_print(10,0,"|",@points);
		&get_kpath_kpf_print(11,7,"-",@points);
		&get_kpath_kpf_print(12,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(13,6,"X",@points);
		} else {
			&get_kpath_kpf_print(13,6,"|",@points);
			&get_kpath_kpf_print(14,0,"-",@points);
			&get_kpath_kpf_print(15,23,"|",@points);
			&get_kpath_kpf_print(16,24,"-",@points);
			&get_kpath_kpf_print(17,16,"-",@points);
			&get_kpath_kpf_print(18,0,"-",@points);
			&get_kpath_kpf_print(19,19,"-",@points);
			&get_kpath_kpf_print(20,26,"|",@points);
			&get_kpath_kpf_print(21,27,"-",@points);
			&get_kpath_kpf_print(22,18,"-",@points);
			&get_kpath_kpf_print(23,0,"|",@points);
			&get_kpath_kpf_print(24,22,"-",@points);
			&get_kpath_kpf_print(25,0,"-",@points);
			&get_kpath_kpf_print(26,21,"X",@points);
		}
	} elsif ($latticetype eq "mC2"){
#mC2:b>asin(beta) BZ 12-face 0-1-3-2-0|5-0-4
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,2,"-",@points);
		&get_kpath_kpf_print(5,0,"|",@points);
		&get_kpath_kpf_print(6,5,"-",@points);
		&get_kpath_kpf_print(7,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(8,4,"X",@points);
		} else {
			&get_kpath_kpf_print(8,4,"|",@points);
			&get_kpath_kpf_print(9,0,"-",@points);
			&get_kpath_kpf_print(10,16,"-",@points);
			&get_kpath_kpf_print(11,18,"-",@points);
			&get_kpath_kpf_print(12,17,"-",@points);
			&get_kpath_kpf_print(13,0,"|",@points);
			&get_kpath_kpf_print(14,20,"-",@points);
			&get_kpath_kpf_print(15,0,"-",@points);
			&get_kpath_kpf_print(16,19,"X",@points);
		}
	} elsif ($latticetype eq "mC3"){
#mC3:b>asin(beta) BZ 14-face 0-2-8-7-3-0-1|6-0-5
		&get_kpath_kpf_print(2,2,"-",@points);
		&get_kpath_kpf_print(3,8,"|",@points);
		&get_kpath_kpf_print(4,7,"-",@points);
		&get_kpath_kpf_print(5,3,"-",@points);
		&get_kpath_kpf_print(6,0,"-",@points);
		&get_kpath_kpf_print(7,1,"|",@points);
		&get_kpath_kpf_print(8,6,"-",@points);
		&get_kpath_kpf_print(9,0,"-",@points);
		if ($inversion eq "Y"){
			&get_kpath_kpf_print(10,5,"X",@points);
		} else {
			&get_kpath_kpf_print(10,5,"|",@points);
			&get_kpath_kpf_print(11,0,"-",@points);
			&get_kpath_kpf_print(12,20,"-",@points);
			&get_kpath_kpf_print(13,26,"|",@points);
			&get_kpath_kpf_print(14,25,"-",@points);
			&get_kpath_kpf_print(15,21,"-",@points);
			&get_kpath_kpf_print(16,0,"-",@points);
			&get_kpath_kpf_print(17,19,"|",@points);
			&get_kpath_kpf_print(18,24,"-",@points);
			&get_kpath_kpf_print(19,0,"-",@points);
			&get_kpath_kpf_print(20,23,"X",@points);
		}
	} elsif ($bravais eq "aP"){
#triclinic obtuse 0-3|2-0-1|7-0-6|5-0-4
		if ($points[15] eq "X"){
			&get_kpath_kpf_print(2,3,"|",@points);
			&get_kpath_kpf_print(3,2,"-",@points);
			&get_kpath_kpf_print(4,0,"-",@points);
			&get_kpath_kpf_print(5,1,"|",@points);
			&get_kpath_kpf_print(6,7,"-",@points);
			&get_kpath_kpf_print(7,0,"-",@points);
			&get_kpath_kpf_print(8,6,"|",@points);
			&get_kpath_kpf_print(9,5,"-",@points);
			&get_kpath_kpf_print(10,0,"-",@points);
			if ($inversion eq "Y"){
				&get_kpath_kpf_print(11,4,"X",@points);
			} else {
				&get_kpath_kpf_print(11,4,"|",@points);
				&get_kpath_kpf_print(12,0,"-",@points);
				&get_kpath_kpf_print(13,10,"|",@points);
				&get_kpath_kpf_print(14,9,"-",@points);
				&get_kpath_kpf_print(15,0,"-",@points);
				&get_kpath_kpf_print(16,8,"|",@points);
				&get_kpath_kpf_print(17,14,"-",@points);
				&get_kpath_kpf_print(18,0,"-",@points);
				&get_kpath_kpf_print(19,13,"|",@points);
				&get_kpath_kpf_print(20,12,"-",@points);
				&get_kpath_kpf_print(21,0,"-",@points);
				&get_kpath_kpf_print(22,11,"X",@points);
			}
		} elsif ($points[15] eq "Y2"){
#triclinic acute 0-4|2-0-1|8-0-7|6-0-5
			&get_kpath_kpf_print(2,4,"|",@points);
			&get_kpath_kpf_print(3,2,"-",@points);
			&get_kpath_kpf_print(4,0,"-",@points);
			&get_kpath_kpf_print(5,1,"|",@points);
			&get_kpath_kpf_print(6,8,"-",@points);
			&get_kpath_kpf_print(7,0,"-",@points);
			&get_kpath_kpf_print(8,7,"|",@points);
			&get_kpath_kpf_print(9,6,"-",@points);
			&get_kpath_kpf_print(10,0,"-",@points);
			if ($inversion eq "Y"){
				&get_kpath_kpf_print(11,5,"X",@points);
			} else {
				&get_kpath_kpf_print(11,5,"|",@points);
				&get_kpath_kpf_print(12,0,"-",@points);
				&get_kpath_kpf_print(13,12,"|",@points);
				&get_kpath_kpf_print(14,10,"-",@points);
				&get_kpath_kpf_print(15,0,"-",@points);
				&get_kpath_kpf_print(16,9,"|",@points);
				&get_kpath_kpf_print(17,16,"-",@points);
				&get_kpath_kpf_print(18,0,"-",@points);
				&get_kpath_kpf_print(19,15,"|",@points);
				&get_kpath_kpf_print(20,14,"-",@points);
				&get_kpath_kpf_print(21,0,"-",@points);
				&get_kpath_kpf_print(22,13,"X",@points);
			}
		} else {
			die ("get_kpath_kpf: triclinic type ga okasii \n");
		}
	} else {
		die ("get_kpath_kpf: extended Bravais ga okasii \n");
	}
	close OUT;
}

# SC VERSION!
#lattice_symbol (CUB nado) wo return de kaesu

#sosite conventional or primitive wo kaesu
sub get_standard_SC{
#kuukangun
	my @t=&get_bposcar;
	my $center=$t[2];
	my $lattice_type=$t[3];
	my $bravais=$t[4];
#lattice symbol
#lattice parameters real
	my @abcreal=&get_abc(@latvec);
#new primitive lattice vector
	my @latvec_new;
	my $lattice_symbol="BUG";
#jiku wo torinaosu, cartesian ni suru
	&d2c; 
	if ($bravais eq "cP"){
		$lattice_symbol="CUB";
		&get_standard_SC_get_latvec("P");
	}
	if ($bravais eq "cF"){
		$lattice_symbol="FCC";
		if ($_[0] eq "primitive"){
			&get_standard_SC_get_latvec("F");
		} else {
			&get_standard_SC_get_latvec("P");
		}
	}
	if ($bravais eq "cI"){
		$lattice_symbol="BCC";
		if ($_[0] eq "primitive"){
			&get_standard_SC_get_latvec("I");
		} else {
			&get_standard_SC_get_latvec("P");
		}
	}
	if ($bravais eq "tP"){
		$lattice_symbol="TET";
		&get_standard_SC_get_latvec("P");
	}
	if ($bravais eq "tI"){
		if ($abcreal[2] < $abcreal[0]){
			$lattice_symbol="BCT1";
		} else {
			$lattice_symbol="BCT2";
		}
		if ($_[0] eq "primitive"){
			&get_standard_SC_get_latvec("I");
		} else {
			&get_standard_SC_get_latvec("P");
		}
	}
	if ($lattice_type eq "Orthorhombic"){
#c-center
		if ($center eq "A"){
			$lattice_symbol="ORCC";
			if ($abcreal[1] < $abcreal[2]) {
				&rotate_axis("bca");
			} else {
				&rotate_axis("cb-a");
			}
			if ($_[0] eq "primitive"){
				&get_standard_SC_get_latvec("OC");
			} else {
				&get_standard_SC_get_latvec("P");
			}
		} elsif ($center eq "C"){
			$lattice_symbol="ORCC";
			if ($abcreal[1] < $abcreal[0]) {
				&rotate_axis("ba-c");
			}
			if ($_[0] eq "primitive"){
				&get_standard_SC_get_latvec("OC");
			} else {
				&get_standard_SC_get_latvec("P");
			}
		} else {
#not c-center: subete onaji kaiten gyouretu ni naru
			if (($abcreal[0] <= $abcreal[1]) && ($abcreal[1] <= $abcreal[2])){
				&rotate_axis("abc");
			} elsif (($abcreal[0] <= $abcreal[2]) && ($abcreal[2] < $abcreal[1])){
				&rotate_axis("-acb");
			} elsif (($abcreal[1] < $abcreal[0]) && ($abcreal[0] <= $abcreal[2])){
				&rotate_axis("ba-c");
			} elsif (($abcreal[1] <= $abcreal[2]) && ($abcreal[2] < $abcreal[0])){
				&rotate_axis("bca");
			} elsif (($abcreal[2] < $abcreal[0]) && ($abcreal[0] <= $abcreal[1])){
				&rotate_axis("cab");
			} elsif (($abcreal[2] < $abcreal[1]) && ($abcreal[1] < $abcreal[0])){
				&rotate_axis("c-ba");
			} else {
				die ("ORC, ORCF, ORCI no izureka, kaiseki de bug: kousiteisuu ga doreka onaji?");
			}
			if ($center eq "P"){
				$lattice_symbol="ORC";
				&get_standard_SC_get_latvec("P");
			}
			if ($center eq "F"){
				@abcreal=&get_abc(@latvec);
				my $tx=1/$abcreal[0]/$abcreal[0];
				my $ty=1/$abcreal[1]/$abcreal[1]+1/$abcreal[2]/$abcreal[2];
				if (abs($tx-$ty) < 0.00001){
					$lattice_symbol="ORCF3";
				} elsif ($tx>$ty){
					$lattice_symbol="ORCF1";
				} else {
					$lattice_symbol="ORCF2";
				}
				if ($_[0] eq "primitive"){
					&get_standard_SC_get_latvec("F");
				} else {
					&get_standard_SC_get_latvec("P");
				}
			}
			if ($center eq "I"){
				$lattice_symbol="ORCI";
				if ($_[0] eq "primitive"){
					&get_standard_SC_get_latvec("I");
				} else {
					&get_standard_SC_get_latvec("P");
				}
			}
		}
	}
	if ($lattice_type eq "Hexagonal_Rhombohedral"){
		if ($center eq "P"){
			$lattice_symbol="HEX";
			&get_standard_SC_get_latvec("H");
		}
		if ($center eq "R"){
#primitive? hex? setting wo kakunin
			my $delta_length=&max($abcreal[0],$abcreal[1],$abcreal[2])-&min($abcreal[0],$abcreal[1],$abcreal[2]);
			my $delta_angle=&max($abcreal[3],$abcreal[4],$abcreal[5])-&min($abcreal[3],$abcreal[4],$abcreal[5]);
			if (($delta_length < 0.001) && ($delta_angle < 0.001)){
#primitive: nani mo shinai
			} elsif ((abs($abcreal[3]-90) < 0.001) && (abs($abcreal[5]-120) < 0.001)){
#hexagonal setting = obverse ni kaeru
				my @a=&sum_vv(&product_vs(@latvec[0..2],2/3),&product_vs(@latvec[3..5],1/3));
				@a=&sum_vv(@a,&product_vs(@latvec[6..8],1/3));
				my @b=&sum_vv(&product_vs(@latvec[0..2],-1/3),&product_vs(@latvec[3..5],1/3));
				@b=&sum_vv(@b,&product_vs(@latvec[6..8],1/3));
				my @c=&sum_vv(&product_vs(@latvec[0..2],-1/3),&product_vs(@latvec[3..5],-2/3));
				@c=&sum_vv(@c,&product_vs(@latvec[6..8],1/3));
				@latvec=(@a,@b,@c);
				@abcreal=&get_abc(@latvec);
			} else {
				die ("Rhombohedral, POSCAR no lattice parameter ga okasii \n");
			} 
			if ($abcreal[3] < 90){
				$lattice_symbol="RHL1";
			} else {
				$lattice_symbol="RHL2";
			}
			get_standard_SC_get_latvec("R");
		}
	}
	if ($lattice_type eq "Monoclinic"){
		if ($center eq "P"){
			if ($abcreal[0] < $abcreal[2]) {
				&rotate_axis("b-ac");
			} else {
				&rotate_axis("-bc-a");
			}
			$lattice_symbol="MCL";
			@latvec_new=&get_standard_SC_get_latvec("P");
		}
		if ($center eq "C"){
		&rotate_axis("b-ac");
#MCLC: need lattice parameters of PRIMITIVE reci: 
			@abcreal=&get_abc(@latvec);
			my $ta=$abcreal[0]/2;
			my $tb=$abcreal[1]/2;
			my $tp=$abcreal[0]*$abcreal[0];
			my $tq=$abcreal[1]*$abcreal[1]*&sind($abcreal[3])*&sind($abcreal[3]);
			my $coskga=($tp-$tq)/($tp+$tq);
#cos 91=-0.0174
			if (abs($coskga) < 0.000175){
				$lattice_symbol="MCLC2";
			} elsif ($coskga < 0){
				$lattice_symbol="MCLC1";
			} else {
				my $tx=$abcreal[1]*&cosd($abcreal[3])/$abcreal[2]+$abcreal[1]*$abcreal[1]*&sind($abcreal[3])*&sind($abcreal[3])/$abcreal[0]/$abcreal[0];
				if (abs($tx-1) < 0.0001){
					$lattice_symbol="MCLC4";
				} elsif ($tx < 1){
					$lattice_symbol="MCLC3";
				} else {
					$lattice_symbol="MCLC5";
				}
			}
			if ($_[0] eq "primitive"){
				&get_standard_SC_get_latvec("MC");
			} else {
				&get_standard_SC_get_latvec("P");
			}
		}
	}
	if ($bravais eq "aP"){
		&c2d;
# get niggli-reci poscar
		my @latvec_reci=&tmatrix3(&invmatrix3(@latvec));
		my @latvec_reci_niggli=&get_niggli_latvec(@latvec_reci);
		my @latvec_niggli=&tmatrix3(&invmatrix3(@latvec_reci_niggli));
		my @conversion_matrix=&round_array(&product_mm(@latvec,&invmatrix3(@latvec_niggli)));
		for (my $i=0; $i<$num_atoms; $i++){
			($x[$i],$y[$i],$z[$i])=&product_vm($x[$i],$y[$i],$z[$i], @conversion_matrix);		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
		}
		@latvec=@latvec_niggli;
		&d2c;
#k_gamma wo 90 do ni chikazukeru
		@abcreal=&get_abc(@latvec);
		my @recilatvec=&reci_latvec(@latvec);
		my @abcreci=&get_abc(@recilatvec);
		@abcreci=&sum_vv(@abcreci,1E-12,2E-12,3E-12,1E-12,2E-12,3E-12);
		my $diffx=abs($abcreci[3]-90);
		my $diffy=abs($abcreci[4]-90);
		my $diffz=abs($abcreci[5]-90);
		if ($diffx == &min($diffx,$diffy,$diffz)){
			&rotate_axis("bca");
		} elsif ($diffy == &min($diffx,$diffy,$diffz)){
			&rotate_axis("cab");
		} 
		@abcreal=&get_abc(@latvec);
		@recilatvec=&reci_latvec(@latvec);
		@abcreci=&get_abc(@recilatvec);
#zenbu >90 matawa <90 ni suru
		if (($abcreci[3] <= 90) && ($abcreci[4] > 90) && ($abcreci[5] > 90)){
			&rotate_axis("a-b-c");
		} elsif (($abcreci[3] > 90) && ($abcreci[4] <= 90) && ($abcreci[5] > 90)){
			&rotate_axis("-ab-c");
		} elsif (($abcreci[3] > 90) && ($abcreci[4] > 90) && ($abcreci[5] <= 90)){
			&rotate_axis("-a-bc");
		} elsif (($abcreci[3] > 90) && ($abcreci[4] <= 90) && ($abcreci[5] <= 90)){
			&rotate_axis("a-b-c");
		} elsif (($abcreci[3] <= 90) && ($abcreci[4] > 90) && ($abcreci[5] <= 90)){
			&rotate_axis("-ab-c");
		} elsif (($abcreci[3] <= 90) && ($abcreci[4] <= 90) && ($abcreci[5] > 90)){
			&rotate_axis("-a-bc");
		}
		@abcreal=&get_abc(@latvec);
		@recilatvec=&reci_latvec(@latvec);
		@abcreci=&get_abc(@recilatvec);
		if ($abcreci[5] > 90.001){
			$lattice_symbol="TRI1a";
		} elsif ($abcreci[5] < 89.999){
			$lattice_symbol="TRI1b";
		} elsif (($abcreci[3] > 90) || ($abcreci[4] > 90)){
			$lattice_symbol="TRI2a";
		} else {
			$lattice_symbol="TRI2b";
		}
		&get_standard_SC_get_latvec("P");
	}
	for (my $i=0; $i<$num_atoms; $i++){
		($x[$i],$y[$i],$z[$i])=&site_in_cell($x[$i],$y[$i],$z[$i]);
	}
	return ($lattice_symbol);
}

# SC VERSION!
# kousi teisuu wo kaku
# P, C, I, F, H, R no latvec wo kaesu
# centering wo hikisuu to site ataeru
# Hexagonal = H not P!
# @latvec wo tsukau
# atarasii latvec wo kaesu
# zahyo ha cartesian to suru! direct de kaesu
sub get_standard_SC_get_latvec{
	my @latvec_new;
	my @abcreal=&get_abc(@latvec);
	my $ta=$abcreal[0]/2;
	my $tb=$abcreal[1]/2;
	my $tc=$abcreal[2]/2;
	if ($_[0] eq "P"){
		&c2d;
		my @t=&clat(@abcreal,"s");
		@latvec_new=($t[0],0,0,$t[1],$t[2],0,@t[3..5]);
	} elsif ($_[0] eq "MC"){
		my @a=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],0.5));
		my @b=&sum_vv(&product_vs(@latvec[0..2],-0.5),&product_vs(@latvec[3..5],0.5));
		@latvec=(@a,@b,@latvec[6..8]);
		&c2d;
		@latvec_new=($ta,$tb,0,-$ta,$tb,0,0,$abcreal[2]*&cosd($abcreal[3]),$abcreal[2]*&sind($abcreal[3]));
	} elsif ($_[0] eq "OC"){
		my @a=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],-0.5));
		my @b=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],0.5));
		@latvec=(@a,@b,@latvec[6..8]);
		&c2d;
		@latvec_new=($ta,-$tb,0,$ta,$tb,0,0,0,$abcreal[2]);
	} elsif ($_[0] eq "I"){
		my @a=&sum_vv(&product_vs(@latvec[0..2],-0.5),&product_vs(@latvec[3..5],0.5));
		@a=&sum_vv(@a,&product_vs(@latvec[6..8],0.5));
		my @b=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],-0.5));
		@b=&sum_vv(@b,&product_vs(@latvec[6..8],0.5));
		my @c=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],0.5));
		@c=&sum_vv(@c,&product_vs(@latvec[6..8],-0.5));
		@latvec=(@a,@b,@c);
		&c2d;
		@latvec_new=(-$ta,$tb,$tc,$ta,-$tb,$tc,$ta,$tb,-$tc);
	} elsif ($_[0] eq "F"){
		my @a=&sum_vv(&product_vs(@latvec[3..5],0.5),&product_vs(@latvec[6..8],0.5));
		my @b=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[6..8],0.5));
		my @c=&sum_vv(&product_vs(@latvec[0..2],0.5),&product_vs(@latvec[3..5],0.5));
		@latvec=(@a,@b,@c);
		&c2d;
		@latvec_new=(0,$tb,$tc,$ta,0,$tc,$ta,$tb,0);
	} elsif ($_[0] eq "H"){
		my $tx=$abcreal[0]*sqrt(3)/2;
		&c2d;
		@latvec_new=($ta,-$tx,0,$ta,$tx,0,0,0,$abcreal[2]);
	} elsif ($_[0] eq "R"){
		my $tx=$abcreal[0]*&cosd($abcreal[3]/2);
		my $ty=$abcreal[0]*&sind($abcreal[3]/2);
		my $tz=&cosd($abcreal[3])/&cosd($abcreal[3]/2);
		&c2d;
		@latvec_new=($tx,-$ty,0,$tx,$ty,0,$abcreal[0]*$tz,0,$abcreal[0]*sqrt(1-$tz*$tz));
	} else {
		die ("Bug in get_standard_SC_get_latvec\n");
	}
	for (my $i=0; $i<=8;$i++){
		$latvec[$i]=&precise($latvec_new[$i],0);
	}
	&unique();
}





#Setyawan Curtarolo ni motozuita
#Brillouin zone no ten
#Shuturyoku ha array ($latticetype, @zahyou, $symbol ...)
#@zahyou: eg (0,0,0), $symbol: eg Gamma, A,
#Kyousei teki ni standard primitive ni suru

sub bz_point_SC{
	my $latticetype=&get_standard_SC("conventional");
	my @abcreal=&get_abc(@latvec);
	my ($a, $b, $c)=@abcreal[0..2];
	my $cal=&cosd($abcreal[3]);
	my $cbe=&cosd($abcreal[4]);
	my $cga=&cosd($abcreal[5]);
	my $sal=&sind($abcreal[3]);
	my $sbe=&sind($abcreal[4]);
	my $sga=&sind($abcreal[5]);
	my @t=($latticetype,0,0,0,"Gamma");
	if ($latticetype eq "CUB"){
		push @t, (1/2,1/2,0,"M");
		push @t, (1/2,1/2,1/2,"R");
		push @t, (0,1/2,0,"X");
	} elsif ($latticetype eq "FCC"){
		push @t, (3/8,3/8,3/4,"K");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (5/8,1/4,5/8,"U");
		push @t, (1/2,1/4,3/4,"W");
		push @t, (1/2,0,1/2,"X");
	} elsif ($latticetype eq "BCC"){
		push @t, (1/2,-1/2,1/2,"H");
		push @t, (1/4,1/4,1/4,"P");
		push @t, (0,0,1/2,"N");
	} elsif ($latticetype eq "TET"){
		push @t, (1/2,1/2,1/2,"A");
		push @t, (1/2,1/2,0,"M");
		push @t, (0,1/2,1/2,"R");
		push @t, (0,1/2,0,"X");
		push @t, (0,0,1/2,"Z");
	} elsif ($latticetype eq "BCT1"){
		my $y=(1+$c*$c/$a/$a)/4;
		push @t, (-1/2,1/2,1/2,"M");
		push @t, (0,1/2,0,"N");
		push @t, (1/4,1/4,1/4,"P");
		push @t, (0,0,1/2,"X");
		push @t, ($y,$y,-$y,"Z");
		push @t, (-$y,1-$y,$y,"Z1");
	} elsif ($latticetype eq "BCT2"){
		my $y=(1+$a*$a/$c/$c)/4;
		my $z=$a*$a/2/$c/$c;
		push @t, (0,1/2,0,"N");
		push @t, (1/4,1/4,1/4,"P");
		push @t, (-$y,$y,$y,"Sigma");
		push @t, ($y,1-$y,-$y,"Sigma1");
		push @t, (0,0,1/2,"X");
		push @t, (-$z,$z,1/2,"Y");
		push @t, (1/2,1/2,-$z,"Y1");
		push @t, (1/2,1/2,-1/2,"Z");
	} elsif ($latticetype eq "ORC"){
		push @t, (1/2,1/2,1/2,"R");
		push @t, (1/2,1/2,0,"S");
		push @t, (0,1/2,1/2,"T");
		push @t, (1/2,0,1/2,"U");
		push @t, (1/2,0,0,"X");
		push @t, (0,1/2,0,"Y");
		push @t, (0,0,1/2,"Z");	
	} elsif (($latticetype eq "ORCF1")||($latticetype eq "ORCF3")){
		my $z=(1+$a*$a/$b/$b-$a*$a/$c/$c)/4;
		my $y=(1+$a*$a/$b/$b+$a*$a/$c/$c)/4;
		push @t, (1/2,1/2+$z,$z,"A");
		push @t, (1/2,1/2-$z,1-$z,"A1");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (1,1/2,1/2,"T");
		push @t, (0,$y,$y,"X");
		push @t, (1,1-$y,1-$y,"X1");
		push @t, (1/2,0,1/2,"Y");
		push @t, (1/2,1/2,0,"Z");
	} elsif ($latticetype eq "ORCF2"){
		my $y=(1+$a*$a/$b/$b-$a*$a/$c/$c)/4;
		my $p=(1+$c*$c/$b/$b-$c*$c/$a/$a)/4;
		my $d=(1+$b*$b/$a/$a-$b*$b/$c/$c)/4;
		push @t, (1/2,1/2-$y,1-$y,"C");
		push @t, (1/2,1/2+$y,$y,"C1");
		push @t, (1/2-$d,1/2,1-$d,"D");
		push @t, (1/2+$d,1/2,$d,"D1");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (1-$p,1/2-$p,1/2,"H");
		push @t, ($p,1/2+$p,1/2,"H1");
		push @t, (0,1/2,1/2,"X");
		push @t, (1/2,0,1/2,"Y");
		push @t, (1/2,1/2,0,"Z");
	} elsif ($latticetype eq "ORCI"){
		my $z=(1+$a*$a/$c/$c)/4;
		my $y=(1+$b*$b/$c/$c)/4;
		my $d=($b*$b-$a*$a)/$c/$c/4;
		my $m=($b*$b+$a*$a)/$c/$c/4;
		push @t, (-$m,$m,1/2-$d,"L");
		push @t, ($m,-$m,1/2+$d,"L1");
		push @t, (1/2-$d,1/2+$d,-$m,"L2");
		push @t, (0,1/2,0,"R");
		push @t, (1/2,0,0,"S");
		push @t, (0,0,1/2,"T");
		push @t, (1/4,1/4,1/4,"W");
		push @t, (-$z,$z,$z,"X");
		push @t, ($z,1-$z,-$z,"X1");
		push @t, ($y,-$y,$y,"Y");
		push @t, (1-$y,$y,-$y,"Y1");
		push @t, (1/2,1/2,-1/2,"Z");
	} elsif ($latticetype eq "ORCC"){
		my $z=(1+$a*$a/$b/$b)/4;
		push @t, ($z,$z,1/2,"A");
		push @t, (-$z,1-$z,1/2,"A1");
		push @t, (0,1/2,13/2,"R");
		push @t, (0,1/2,0,"S");
		push @t, (-1/2,1/2,1/2,"T");
		push @t, ($z,$z,0,"X");
		push @t, (-$z,1-$z,0,"X1");
		push @t, (-1/2,1/2,0,"Y");
		push @t, (0,0,1/2,"Z");
	} elsif ($latticetype eq "HEX"){
		push @t, (0,0,1/2,"A");
		push @t, (1/3,1/3,1/2,"H");
		push @t, (1/3,1/3,0,"K");
		push @t, (1/2,0,1/2,"L");
		push @t, (1/2,0,0,"M");
	} elsif ($latticetype eq "RHL1"){
		my $y=(1+4*$cal)/(2+4*$cal);
		my $n=3/4-$y/2;
		push @t, ($y,1/2,1-$y,"B");
		push @t, (1/2,1-$y,$y-1,"B1");
		push @t, (1/2,1/2,0,"F");
		push @t, (1/2,0,0,"L");
		push @t, (0,0,-1/2,"L1");
		push @t, ($y,$n,$n,"P");
		push @t, (1-$n,1-$n,1-$y,"P1");
		push @t, ($n,$n,$y-1,"P2");
		push @t, (1-$n,$n,0,"Q");
		push @t, ($n,0,-$n,"X");
		push @t, (1/2,1/2,1/2,"Z");
	} elsif ($latticetype eq "RHL2"){
		my $y=1/2/&tand($abcreal[3]/2)/&tand($abcreal[3]/2);
		my $n=3/4-$y/2;
		push @t, (1/2,-1/2,0,"F");
		push @t, (1/2,0,0,"L");
		push @t, (1-$n,-$n,1-$n,"P");
		push @t, ($n,$n-1,$n-1,"P1");
		push @t, ($y,$y,$y,"Q");
		push @t, (1-$y,-$y,-$y,"Q1");
		push @t, (1/2,-1/2,1/2,"Z");
	} elsif ($latticetype eq "MCL"){
		my $y=(1-$b*$cal/$c)/(2*$sal*$sal);
		my $n=1/2-$y*$c*$cal/$b;
		push @t, (1/2,1/2,0,"A");
		push @t, (0,1/2,1/2,"C");
		push @t, (1/2,0,1/2,"D");
		push @t, (1/2,0,-1/2,"D1");
		push @t, (1/2,1/2,1/2,"E");
		push @t, (0,$y,1-$n,"H");
		push @t, (0,1-$y,$n,"H1");
		push @t, (0,$y,-$n,"H2");
		push @t, (1/2,$y,1-$n,"M");
		push @t, (1/2,1-$y,$n,"M1");
		push @t, (1/2,$y,-$n,"M2");
		push @t, (0,1/2,0,"X");
		push @t, (0,0,1/2,"Y");			
		push @t, (0,0,-1/2,"Y1");			
		push @t, (1/2,0,0,"Z");
	} elsif (($latticetype eq "MCLC1")||($latticetype eq "MCLC2")){
		my $z=(2-$b*$cal/$c)/4/$sal/$sal;
		my $y=1/2+2*$z*$c*$cal/$b;
		my $s=3/4-$a*$a/4/$b/$b/$sal/$sal;
		my $p=$s+(3/4-$s)*$b*$cal/$c;
		push @t, (1/2,0,0,"N");
		push @t, (0,-1/2,0,"N1");
		push @t, (1-$z,1-$z,1-$y,"F");
		push @t, ($z,$z,$y,"F1");
		push @t, (-$z,-$z,1-$y,"F2");
		push @t, (1-$z,-$z,1-$y,"F3");
		push @t, ($p,1-$p,1/2,"I");
		push @t, (1-$p,$p-1,1/2,"I1");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (1/2,0,1/2,"M");
		push @t, (1-$s,$s-1,0,"X");
		push @t, ($s,1-$s,0,"X1");
		push @t, ($s-1,-$s,0,"X2");
		push @t, (1/2,1/2,0,"Y");
		push @t, (-1/2,-1/2,0,"Y1");
		push @t, (0,0,1/2,"Z");
	} elsif (($latticetype eq "MCLC3")||($latticetype eq "MCLC4")){
		my $m=(1+$b*$b/$a/$a)/4;
		my $d=$b*$c*$cal/2/$a/$a;
		my $z=$m-1/4+(1-$b*$cal/$c)/4/$sal/$sal;
		my $y=1/2+2*$z*$c*$cal/$b;
		my $p=1+$z-2*$m;
		my $s=$y-2*$d;
		push @t, (1-$p,1-$p,1-$s,"F");
		push @t, ($p,$p-1,$s,"F1");
		push @t, (1-$p,-$p,1-$s,"F2");
		push @t, ($z,$z,$y,"H");
		push @t, (1-$z,-$z,1-$y,"H1");
		push @t, (-$z,-$z,1-$y,"H2");
		push @t, (1/2,-1/2,1/2,"I");
		push @t, (1/2,0,1/2,"M");
		push @t, (1/2,0,0,"N");
		push @t, (0,-1/2,0,"N1");
		push @t, (1/2,-1/2,0,"X");
		push @t, ($m,$m,$d,"Y");
		push @t, (1-$m,-$m,-$d,"Y1");
		push @t, (-$m,-$m,-$d,"Y2");
		push @t, ($m,$m-1,$d,"Y3");
		push @t, (0,0,1/2,"Z");
	} elsif ($latticetype eq "MCLC5"){
		my $z=($b*$b/$a/$a+(1-$b*$cal/$c)/$sal/$sal)/4;
		my $y=1/2+2*$z*$c*$cal/$b;
		my $m=$y/2+$b*$b/4/$a/$a-$b*$c*$cal/2/$a/$a;
		my $n=2*$m-$z;
		my $o=(4*$n-1-$b*$b*$sal*$sal/$a/$a)*$c/(2*$b*$cal);
		my $d=$z*$c*$cal/$b+$o/2-1/4;
		my $r=1-$z*$a*$a/$b/$b;
		push @t, ($n,$n,$o,"F");
		push @t, (1-$n,1-$n,1-$o,"F1");
		push @t, ($n,$n-1,$o,"F2");
		push @t, ($z,$z,$y,"H");
		push @t, (1-$z,-$z,1-$y,"H1");
		push @t, (-$z,-$z,1-$y,"H2");
		push @t, ($r,1-$r,1/2,"I");
		push @t, (1-$r,$r-1,1/2,"I1");
		push @t, (1/2,1/2,1/2,"L");
		push @t, (1/2,0,1/2,"M");
		push @t, (1/2,0,0,"N");
		push @t, (0,-1/2,0,"N1");
		push @t, (1/2,-1/2,0,"X");
		push @t, ($m,$m,$d,"Y");
		push @t, (1-$m,-$m,-$d,"Y1");
		push @t, (-$m,-$m,-$d,"Y2");
		push @t, ($m,$m-1,$d,"Y3");
		push @t, (0,0,1/2,"Z");
	} elsif (($latticetype eq "TRI1a")||($latticetype eq "TRI2a")){
		push @t, (1/2,1/2,0,"L");
		push @t, (0,1/2,1/2,"M");
		push @t, (1/2,0,1/2,"N");
		push @t, (1/2,1/2,1/2,"R");
		push @t, (1/2,0,0,"X");
		push @t, (0,1/2,0,"Y");
		push @t, (0,0,1/2,"Z");
	} elsif (($latticetype eq "TRI1b")||($latticetype eq "TRI2b")){
		push @t, (1/2,-1/2,0,"L");
		push @t, (0,0,1/2,"M");
		push @t, (-1/2,-1/2,1/2,"N");
		push @t, (0,-1/2,1/2,"R");
		push @t, (0,-1/2,0,"X");
		push @t, (1/2,0,0,"Y");
		push @t, (-1/2,0,1/2,"Z");
	} else {
		die ("bz_point_SC: lattice type ga okasii \n");
	}
	&get_standard_SC("primitive");
	return (@t);
}

#Setyawan & Curtarolo no kpath wo motomeru
#kpf huu ni tsubo_temp_get_kpath ni syuturyoku
#file mei ha betu no basyo de kaeru!
sub get_kpath_kpf_SC{
#Conventional wo nyuusyu: k-point data wo motomeru
	my @points=&bz_point_SC;
	my $a=shift @points;
	open OUT, ">tsubo_temp_get_kpath";
	print OUT ("Real form of k-point coordinates (kx,ky,kz,label):\n");
	if ($a eq "CUB"){
#CUB: Gamma-X-M-Gamma-R-X|M-R
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,2,"-",@points);
		&get_kpath_kpf_print(6,3,"|",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		&get_kpath_kpf_print(8,2,"X",@points);
	} elsif ($a eq "FCC"){
#FCC: Gamma-X-W-K-Gamma-L-U-W-L-K|U-X
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,4,"-",@points);
		&get_kpath_kpf_print(4,1,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,2,"-",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,4,"-",@points);
		&get_kpath_kpf_print(9,2,"-",@points);
		&get_kpath_kpf_print(10,1,"|",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		&get_kpath_kpf_print(12,5,"X",@points);
	} elsif ($a eq "BCC"){
#BCC: Gamma-H-N-Gamma-P-H|P-N
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,1,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,2,"-",@points);
		&get_kpath_kpf_print(6,1,"|",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		&get_kpath_kpf_print(8,3,"X",@points);
	} elsif ($a eq "TET"){
#TET: Gamma-X-M-Gamma-Z-R-A-Z|X-R|M-A
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,4,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,5,"-",@points);
		&get_kpath_kpf_print(6,3,"-",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		&get_kpath_kpf_print(8,5,"|",@points);
		&get_kpath_kpf_print(9,4,"-",@points);
		&get_kpath_kpf_print(10,3,"|",@points);
		&get_kpath_kpf_print(11,2,"-",@points);
		&get_kpath_kpf_print(12,1,"X",@points);
	} elsif ($a eq "BCT1"){
#BCT1: Gamma-X-M-Gamma-Z-P-N-Z1-M|X-P
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,4,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,5,"-",@points);
		&get_kpath_kpf_print(6,3,"-",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		&get_kpath_kpf_print(8,6,"-",@points);
		&get_kpath_kpf_print(9,1,"|",@points);
		&get_kpath_kpf_print(10,4,"-",@points);
		&get_kpath_kpf_print(11,3,"X",@points);
	} elsif ($a eq "BCT2"){
#BCT2: Gamma-X-Y-Sigma-Gamma-Z-Sigma1-N-P-Y1-Z|X-P
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,6,"-",@points);
		&get_kpath_kpf_print(4,3,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,4,"-",@points);
		&get_kpath_kpf_print(8,1,"-",@points);
		&get_kpath_kpf_print(9,2,"-",@points);
		&get_kpath_kpf_print(10,7,"-",@points);
		&get_kpath_kpf_print(11,8,"|",@points);
		&get_kpath_kpf_print(12,5,"-",@points);
		&get_kpath_kpf_print(13,2,"X",@points);
	} elsif ($a eq "ORC"){
#ORC: Gamma-X-S-Y-Gamma-Z-U-R-T-Z|Y-T|U-X|S-R
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,2,"-",@points);
		&get_kpath_kpf_print(4,6,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,7,"-",@points);
		&get_kpath_kpf_print(7,4,"-",@points);
		&get_kpath_kpf_print(8,1,"-",@points);
		&get_kpath_kpf_print(9,3,"-",@points);
		&get_kpath_kpf_print(10,7,"|",@points);
		&get_kpath_kpf_print(11,6,"-",@points);
		&get_kpath_kpf_print(12,3,"|",@points);
		&get_kpath_kpf_print(13,4,"-",@points);
		&get_kpath_kpf_print(14,5,"|",@points);
		&get_kpath_kpf_print(15,2,"-",@points);
		&get_kpath_kpf_print(16,1,"X",@points);
	} elsif ($a eq "ORCF1"){
#ORCF1: Gamma-Y-T-Z-Gamma-X-A1-Y|T-X1|X-A-Z|L-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,7,"-",@points);
		&get_kpath_kpf_print(3,4,"-",@points);
		&get_kpath_kpf_print(4,8,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,5,"-",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		&get_kpath_kpf_print(8,7,"|",@points);
		&get_kpath_kpf_print(9,4,"-",@points);
		&get_kpath_kpf_print(10,6,"|",@points);
		&get_kpath_kpf_print(11,5,"-",@points);
		&get_kpath_kpf_print(12,1,"-",@points);
		&get_kpath_kpf_print(13,8,"|",@points);
		&get_kpath_kpf_print(14,3,"-",@points);
		&get_kpath_kpf_print(15,0,"X",@points);
	} elsif (($a eq "ORCF1")||($a eq "ORCF3")){
#ORCF3: Gamma-Y-T-Z-Gamma-X-A1-Y|X-A-Z|L-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,7,"-",@points);
		&get_kpath_kpf_print(3,4,"-",@points);
		&get_kpath_kpf_print(4,8,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,5,"-",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		&get_kpath_kpf_print(8,7,"|",@points);
		&get_kpath_kpf_print(9,5,"-",@points);
		&get_kpath_kpf_print(10,1,"-",@points);
		&get_kpath_kpf_print(11,8,"|",@points);
		&get_kpath_kpf_print(12,3,"-",@points);
		&get_kpath_kpf_print(13,0,"X",@points);
	} elsif ($a eq "ORCF2"){
#ORCF2: Gamma-Y-C-D-X-Gamma-Z-D1-H-C|C1-Z|X-H1|H-Y|L-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,9,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,3,"-",@points);
		&get_kpath_kpf_print(5,8,"-",@points);
		&get_kpath_kpf_print(6,0,"-",@points);
		&get_kpath_kpf_print(7,10,"-",@points);
		&get_kpath_kpf_print(8,4,"-",@points);
		&get_kpath_kpf_print(9,6,"-",@points);
		&get_kpath_kpf_print(10,1,"|",@points);
		&get_kpath_kpf_print(11,2,"-",@points);
		&get_kpath_kpf_print(12,10,"|",@points);
		&get_kpath_kpf_print(13,8,"-",@points);
		&get_kpath_kpf_print(14,7,"|",@points);
		&get_kpath_kpf_print(15,6,"-",@points);
		&get_kpath_kpf_print(16,9,"|",@points);
		&get_kpath_kpf_print(17,5,"-",@points);
		&get_kpath_kpf_print(18,0,"X",@points);
	} elsif ($a eq "ORCI"){
#ORCI: Gamma-X-L-T-W-R-X1-Z-Gamma-Y-S-W|L1-Y|Y1-Z
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,8,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,6,"-",@points);
		&get_kpath_kpf_print(5,7,"-",@points);
		&get_kpath_kpf_print(6,4,"-",@points);
		&get_kpath_kpf_print(7,9,"-",@points);
		&get_kpath_kpf_print(8,12,"-",@points);
		&get_kpath_kpf_print(9,0,"-",@points);
		&get_kpath_kpf_print(10,10,"-",@points);
		&get_kpath_kpf_print(11,5,"-",@points);
		&get_kpath_kpf_print(12,7,"|",@points);
		&get_kpath_kpf_print(13,2,"-",@points);
		&get_kpath_kpf_print(14,10,"|",@points);
		&get_kpath_kpf_print(15,11,"-",@points);
		&get_kpath_kpf_print(16,12,"X",@points);
	} elsif ($a eq "ORCC"){
#ORCC: Gamma-X-S-R-A-Z-Gamma-Y-X1-A1-T-Y|Z-T
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,6,"-",@points);
		&get_kpath_kpf_print(3,4,"-",@points);
		&get_kpath_kpf_print(4,3,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,9,"-",@points);
		&get_kpath_kpf_print(7,0,"-",@points);
		&get_kpath_kpf_print(8,8,"-",@points);
		&get_kpath_kpf_print(9,7,"-",@points);
		&get_kpath_kpf_print(10,2,"-",@points);
		&get_kpath_kpf_print(11,5,"-",@points);
		&get_kpath_kpf_print(12,8,"|",@points);
		&get_kpath_kpf_print(13,9,"-",@points);
		&get_kpath_kpf_print(14,5,"X",@points);
	} elsif ($a eq "HEX"){
#HEX: Gamma-M-K-Gamma-A-L-H-A|L-M|K-H
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,5,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,0,"-",@points);
		&get_kpath_kpf_print(5,1,"-",@points);
		&get_kpath_kpf_print(6,4,"-",@points);
		&get_kpath_kpf_print(7,2,"-",@points);
		&get_kpath_kpf_print(8,1,"|",@points);
		&get_kpath_kpf_print(9,4,"-",@points);
		&get_kpath_kpf_print(10,5,"|",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		&get_kpath_kpf_print(12,2,"X",@points);
	} elsif ($a eq "RHL1"){
#RHL1: Gamma-L-B1|B-Z-Gamma-X|Q-F-P1-Z|L-P
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,4,"-",@points);
		&get_kpath_kpf_print(3,2,"|",@points);
		&get_kpath_kpf_print(4,1,"-",@points);
		&get_kpath_kpf_print(5,11,"-",@points);
		&get_kpath_kpf_print(6,0,"-",@points);
		&get_kpath_kpf_print(7,10,"|",@points);
		&get_kpath_kpf_print(8,9,"-",@points);
		&get_kpath_kpf_print(9,3,"-",@points);
		&get_kpath_kpf_print(10,7,"-",@points);
		&get_kpath_kpf_print(11,11,"|",@points);
		&get_kpath_kpf_print(12,4,"-",@points);
		&get_kpath_kpf_print(13,6,"X",@points);
	} elsif ($a eq "RHL2"){
#RHL2: Gamma-P-Z-Q-Gamma-F-P1-Q1-L-Z
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,3,"-",@points);
		&get_kpath_kpf_print(3,7,"-",@points);
		&get_kpath_kpf_print(4,5,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,1,"-",@points);
		&get_kpath_kpf_print(7,4,"-",@points);
		&get_kpath_kpf_print(8,6,"-",@points);
		&get_kpath_kpf_print(9,2,"-",@points);
		&get_kpath_kpf_print(10,7,"X",@points);
	} elsif ($a eq "MCL"){
#MCL wikipedia ver
#MCL:Gamma-Y-H-C-E-M1-A-X-Gamma-Z-D-M|Z-A|D-Y|X-H1
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,13,"-",@points);
		&get_kpath_kpf_print(3,6,"-",@points);
		&get_kpath_kpf_print(4,2,"-",@points);
		&get_kpath_kpf_print(5,5,"-",@points);
		&get_kpath_kpf_print(6,10,"-",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		&get_kpath_kpf_print(8,12,"-",@points);
		&get_kpath_kpf_print(9,0,"-",@points);
		&get_kpath_kpf_print(10,15,"-",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		&get_kpath_kpf_print(12,9,"|",@points);
		&get_kpath_kpf_print(13,15,"-",@points);
		&get_kpath_kpf_print(14,1,"|",@points);
		&get_kpath_kpf_print(15,3,"-",@points);
		&get_kpath_kpf_print(16,13,"|",@points);
		&get_kpath_kpf_print(17,12,"-",@points);
		&get_kpath_kpf_print(18,7,"X",@points);
=pod
#ronbun ver
#MCL: Gamma-Y-H-C-E-M1-A-X-H1|M-D-Z|Y-D
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,13,"-",@points);
		&get_kpath_kpf_print(3,6,"-",@points);
		&get_kpath_kpf_print(4,2,"-",@points);
		&get_kpath_kpf_print(5,5,"-",@points);
		&get_kpath_kpf_print(6,10,"-",@points);
		&get_kpath_kpf_print(7,1,"-",@points);
		&get_kpath_kpf_print(8,12,"-",@points);
		&get_kpath_kpf_print(9,7,"|",@points);
		&get_kpath_kpf_print(10,9,"-",@points);
		&get_kpath_kpf_print(11,3,"-",@points);
		&get_kpath_kpf_print(12,15,"|",@points);
		&get_kpath_kpf_print(13,13,"-",@points);
		&get_kpath_kpf_print(14,3,"X",@points);
=cut
	} elsif ($a eq "MCLC1"){
#MCLC1 wikipedia ver
#MCLC1: Gamma-Y-F-L-I|I1-Z-Gamma-X|X1-Y|M-Gamma-N|Z-F1
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,16,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,11,"|",@points);
		&get_kpath_kpf_print(10,12,"-",@points);
		&get_kpath_kpf_print(11,14,"|",@points);
		&get_kpath_kpf_print(12,10,"-",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		&get_kpath_kpf_print(14,1,"|",@points);
		&get_kpath_kpf_print(15,16,"-",@points);
		&get_kpath_kpf_print(16,4,"X",@points);
=pod
#MCLC1 ronbun ver
#MCLC1: Gamma-Y-F-L-I|I1-Z-F1|Y-X1|X-Gamma-N|M-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,16,"-",@points);
		&get_kpath_kpf_print(8,4,"|",@points);
		&get_kpath_kpf_print(9,14,"-",@points);
		&get_kpath_kpf_print(10,12,"|",@points);
		&get_kpath_kpf_print(11,11,"-",@points);
		&get_kpath_kpf_print(12,0,"-",@points);
		&get_kpath_kpf_print(13,1,"|",@points);
		&get_kpath_kpf_print(14,10,"-",@points);
		&get_kpath_kpf_print(15,0,"X",@points);
=cut
	} elsif ($a eq "MCLC2"){
#MCLC2 wikipedia ver
#MCLC2: Gamma-Y-F-L-I|I1-Z-Gamma-M|N-Gamma|Z-F1
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,16,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,10,"|",@points);
		&get_kpath_kpf_print(12,1,"-",@points);
		&get_kpath_kpf_print(13,0,"|",@points);
		&get_kpath_kpf_print(15,16,"-",@points);
		&get_kpath_kpf_print(16,4,"X",@points);
=pod
#MCLC2 ronbun ver
#MCLC2: Gamma-Y-F-L-I|I1-Z-F1|N-Gamma-M
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,3,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,16,"-",@points);
		&get_kpath_kpf_print(8,4,"|",@points);
		&get_kpath_kpf_print(9,1,"-",@points);
		&get_kpath_kpf_print(10,0,"-",@points);
		&get_kpath_kpf_print(11,10,"X",@points);
=cut
	} elsif ($a eq "MCLC3"){
#MCLC3 wikipedia ver
#MCLC3: Gamma-Y-F-H-Z-I-X-Gamma-Z|M-Gamma-N|X-Y1-H1|I-F1
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,12,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,16,"-",@points);
		&get_kpath_kpf_print(6,7,"-",@points);
		&get_kpath_kpf_print(7,11,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,16,"|",@points);
		&get_kpath_kpf_print(10,8,"-",@points);
		&get_kpath_kpf_print(11,0,"-",@points);
		&get_kpath_kpf_print(12,9,"|",@points);
		&get_kpath_kpf_print(13,11,"-",@points);
		&get_kpath_kpf_print(14,13,"-",@points);
		&get_kpath_kpf_print(15,5,"|",@points);
		&get_kpath_kpf_print(16,7,"-",@points);
		&get_kpath_kpf_print(17,2,"X",@points);
=pod
#ronbun ver
#MCLC3: Gamma-Y-F-H-Z-I-F1|H1-Y1-X-Gamma-N|M-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,12,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,16,"-",@points);
		&get_kpath_kpf_print(6,7,"-",@points);
		&get_kpath_kpf_print(7,2,"|",@points);
		&get_kpath_kpf_print(8,5,"-",@points);
		&get_kpath_kpf_print(9,13,"-",@points);
		&get_kpath_kpf_print(10,11,"-",@points);
		&get_kpath_kpf_print(11,0,"-",@points);
		&get_kpath_kpf_print(12,9,"|",@points);
		&get_kpath_kpf_print(13,8,"-",@points);
		&get_kpath_kpf_print(14,0,"X",@points);
=cut
	} elsif ($a eq "MCLC4"){
#MCLC4 wikipedia ver
#MCLC4: Gamma-Y-F-H-Z-I-X-Gamma-Z|M-Gamma-N|X-Y-H1
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,12,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,16,"-",@points);
		&get_kpath_kpf_print(6,7,"-",@points);
		&get_kpath_kpf_print(7,11,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,16,"|",@points);
		&get_kpath_kpf_print(10,8,"-",@points);
		&get_kpath_kpf_print(11,0,"-",@points);
		&get_kpath_kpf_print(12,9,"|",@points);
		&get_kpath_kpf_print(13,11,"-",@points);
		&get_kpath_kpf_print(14,12,"-",@points);
		&get_kpath_kpf_print(15,5,"X",@points);
=pod
#ronbun ver
#MCLC4: Gamma-Y-F-H-Z-I1|H1-Y1-X-Gamma-N|M-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,12,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,4,"-",@points);
		&get_kpath_kpf_print(5,16,"-",@points);
		&get_kpath_kpf_print(6,7,"|",@points);
		&get_kpath_kpf_print(7,5,"-",@points);
		&get_kpath_kpf_print(8,13,"-",@points);
		&get_kpath_kpf_print(9,11,"-",@points);
		&get_kpath_kpf_print(10,0,"-",@points);
		&get_kpath_kpf_print(11,9,"|",@points);
		&get_kpath_kpf_print(12,8,"-",@points);
		&get_kpath_kpf_print(13,0,"X",@points);
=cut
	} elsif ($a eq "MCLC5"){
#MCLC5 wikipedia ver
#MCLC5: Gamma-Y-F-L-I|I1-Z-Gamma-X-Y1-H1|H-F1|F2-X|M-Gamma-N|H-Z
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,18,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,13,"-",@points);
		&get_kpath_kpf_print(10,15,"-",@points);
		&get_kpath_kpf_print(11,5,"|",@points);
		&get_kpath_kpf_print(12,4,"-",@points);
		&get_kpath_kpf_print(13,2,"|",@points);
		&get_kpath_kpf_print(14,3,"-",@points);
		&get_kpath_kpf_print(15,13,"|",@points);
		&get_kpath_kpf_print(16,10,"-",@points);
		&get_kpath_kpf_print(17,0,"-",@points);
		&get_kpath_kpf_print(18,11,"|",@points);
		&get_kpath_kpf_print(19,4,"-",@points);
		&get_kpath_kpf_print(20,18,"X",@points);
=pod
#ronbun ver
#MCLC5: Gamma-Y-F-L-I|I1-Z-H-F1|H1-Y1-X-Gamma-N|M-Gamma
		&get_kpath_kpf_print(1,0,"-",@points);
		&get_kpath_kpf_print(2,14,"-",@points);
		&get_kpath_kpf_print(3,1,"-",@points);
		&get_kpath_kpf_print(4,9,"-",@points);
		&get_kpath_kpf_print(5,7,"|",@points);
		&get_kpath_kpf_print(6,8,"-",@points);
		&get_kpath_kpf_print(7,18,"-",@points);
		&get_kpath_kpf_print(8,4,"-",@points);
		&get_kpath_kpf_print(9,2,"|",@points);
		&get_kpath_kpf_print(10,5,"-",@points);
		&get_kpath_kpf_print(11,15,"-",@points);
		&get_kpath_kpf_print(12,13,"-",@points);
		&get_kpath_kpf_print(13,0,"-",@points);
		&get_kpath_kpf_print(14,11,"|",@points);
		&get_kpath_kpf_print(15,10,"-",@points);
		&get_kpath_kpf_print(16,0,"X",@points);
=cut
	} else {
#Triclinic: X-Gamma-Y|L-Gamma-Z|N-Gamma-M|R-Gamma
		&get_kpath_kpf_print(1,5,"-",@points);
		&get_kpath_kpf_print(2,0,"-",@points);
		&get_kpath_kpf_print(3,6,"|",@points);
		&get_kpath_kpf_print(4,1,"-",@points);
		&get_kpath_kpf_print(5,0,"-",@points);
		&get_kpath_kpf_print(6,7,"|",@points);
		&get_kpath_kpf_print(7,3,"-",@points);
		&get_kpath_kpf_print(8,0,"-",@points);
		&get_kpath_kpf_print(9,2,"|",@points);
		&get_kpath_kpf_print(10,4,"-",@points);
		&get_kpath_kpf_print(11,0,"X",@points);
	}
	close OUT;
}
