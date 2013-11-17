#!/usr/bin/perl

sub readpwm{
    my (@apwm);
    open(FP,$_[0]);
    while($line=<FP>){
        $line=~ s/\s+$//;
        @fields=split(/\t/,$line);
        my $caption=$fields[0];
        $caption =~ s/^\s+//;
        $caption =~ s/\s+$//;
        my $len=$fields[1];
        $len =~ s/^\s+//;
        $len =~ s/\s+$//;
        my (@pwm);
        my (@tmppwm);
        
        for ($i=0;$i<4;$i++){
            $line=<FP>;
            $line=~ s/\s+$//;
            @fields=split(/\t/,$line);
            for ($j=0;$j<=$#fields;$j++){
                 $fields[$j]=~ s/^\s+//;
                 $fields[$j]=~ s/\s+$//;
                 push(@pwm,$fields[$j]);
            }
        }
        push(@tmppwm,$caption);
        push(@tmppwm,$len);
        for ($i=0;$i<$len;$i++){
             push(@tmppwm,$pwm[$i]);
             push(@tmppwm,$pwm[$i+$len]);
             push(@tmppwm,$pwm[$i+2*$len]);
             push(@tmppwm,$pwm[$i+3*$len]);
        }
        
       # $_[1]=$len;
       # $_[2]=$caption;
        push @apwm,[@tmppwm];
        #reverse complement
        undef(@tmppwm);
        push(@tmppwm,$caption."_rv");
        push(@tmppwm,$len);
        for ($i=$len-1;$i>=0;$i--){
             push(@tmppwm,$pwm[$i+$len]);
             push(@tmppwm,$pwm[$i]);
             push(@tmppwm,$pwm[$i+3*$len]);
             push(@tmppwm,$pwm[$i+2*$len]);
        }
        push @apwm,[@tmppwm];
    }
    close(FP);
    return @apwm;
}
sub showpwm{
    my ($pwm)=$_[0];
    my ($sh)=$_[1];
    $len=(($#{$pwm}+1)-2)/4;
    my $i;
    for ($i=0;$i<4;$i++){
         for ($j=0;$j<$len;$j++){
              print $pwm->[$sh+$j*4+$i]."\t";
         }
         print "\n";
    }    
}


@pwm1=readpwm(@ARGV[0]);
@pwm2=readpwm(@ARGV[1]);

@BK=(0.307,0.307,0.193,0.193);
for ($item1=0;$item1<=$#pwm2;$item1++){
     $w=7;
     $ukpwm=$pwm2[$item1];
     $len2=$ukpwm->[1];
     $cap2=$ukpwm->[0];
     $max_motif=-1.0;
     $max_cap="";
     @max1=();
     @max2=();
     $grade=0.0;
     for ($item2=0;$item2<=$#pwm1;$item2++){
          #print $TFalist[$ia];
          #print ">>>>>".$item1."\t".$item2."\n";
          $kpwm=$pwm1[$item2];
          $len1=$kpwm->[1];
          $cap1=$kpwm->[0];
          if ($len1>=$len2) {
              $shortlen=$len2;
          }
          else {
              $shortlen=$len1;
          }
          if ($shortlen<$w){
              $w=$shortlen-2;
              if ($w<=0) {
                  die;
              }
          }
          
          
          $p1=$len2-$w;
          $step=$len1+$len2-2*$w;
        
          #print $len1."\t".$len2."\t".$step."\t".$w."\n";
          @tmppwm1=();
          @tmppwm2=();
          
          $max_sim=-1.0;
          @max_pwm1=();
          @max_pwm2=();
          $max_shift=0;
          for ($i=0;$i<=$step;$i++){
               $p2=$i;
               if ($p2<$p1){
                   for ($j=0;$j<($p1-$p2);$j++){
                        push(@tmppwm1,@BK);
                   }
               }
               elsif ($p2>$p1){
                        for ($j=0;$j<($p2-$p1);$j++){
                             push(@tmppwm2,@BK);
                        }
               }
               @tpwm1=@$kpwm;
             
               @tpwm2=@$ukpwm;
               shift(@tpwm1);
               shift(@tpwm1);
               shift(@tpwm2);
               shift(@tpwm2);
               push(@tmppwm1,@tpwm1);
               push(@tmppwm2,@tpwm2);
               if (($p2+$len2-1)<($p1+$len1-1)){
                   for ($j=0;$j<($p1+$len1-$p2-$len2);$j++){
                        push(@tmppwm2,@BK);
                   }
               }
               elsif (($p2+$len2-1)>($p1+$len1-1)){
                   for ($j=0;$j<($p2+$len2-$p1-$len1);$j++){
                        push(@tmppwm1,@BK);
                   }
               }
               #print $p2."\n";
               $sum=0.0;
               for ($j=0;$j<=$#tmppwm1;$j+=4){
                    $v=0.0;
                    #print "\n-\n";
                    #print $tmppwm1[$j]."-".$tmppwm2[$j]."\t".(($tmppwm1[$j]-$tmppwm2[$j])*($tmppwm1[$j]-$tmppwm2[$j]))."\n";
                    $v+=(($tmppwm1[$j]-$tmppwm2[$j])*($tmppwm1[$j]-$tmppwm2[$j]));
                    #print $tmppwm1[$j+1]."-".$tmppwm2[$j+1]."\t".(($tmppwm1[$j+1]-$tmppwm2[$j+1])*($tmppwm1[$j+1]-$tmppwm2[$j+1]))."\n";
                    $v+=(($tmppwm1[$j+1]-$tmppwm2[$j+1])*($tmppwm1[$j+1]-$tmppwm2[$j+1]));
                    #print $tmppwm1[$j+2]."-".$tmppwm2[$j+2]."\t".(($tmppwm1[$j+2]-$tmppwm2[$j+2])*($tmppwm1[$j+2]-$tmppwm2[$j+2]))."\n";
                    $v+=(($tmppwm1[$j+2]-$tmppwm2[$j+2])*($tmppwm1[$j+2]-$tmppwm2[$j+2]));
                    #print $tmppwm1[$j+3]."-".$tmppwm2[$j+3]."\t".(($tmppwm1[$j+3]-$tmppwm2[$j+3])*($tmppwm1[$j+3]-$tmppwm2[$j+3]))."\n";
                    #print "--\n";
                    $v+=(($tmppwm1[$j+3]-$tmppwm2[$j+3])*($tmppwm1[$j+3]-$tmppwm2[$j+3]));
                    #print sqrt($v)."\t";
          
                    $sum+=sqrt($v);
               }
               $wlen=($#tmppwm1+1)/4;
               #print "w".$wlen."\n";
               $sum/=$wlen;
               $sum/=sqrt(2);
               #print $sum."\n";
               #print ">".$i."\t".$sum."\n";
               #print "pwm1\n";
               #showpwm(\@tmppwm1,0);
               #print "pwm2\n";
               #showpwm(\@tmppwm2,0);
               if ($max_sim<(1.0-$sum)){
                   $max_sim=(1.0-$sum);
                   $max_shift=$i;
                   undef @max_pwm1;
                   undef @max_pwm2;
                   push(@max_pwm1,@tmppwm1);
                   push(@max_pwm2,@tmppwm2);
               }
               undef @tmppwm1;
               undef @tmppwm2;
               undef @tpwm1;
               undef @tpwm2;
               #print "-\n";
          }
          #print $max_shift."\n";
          if ($max_motif<$max_sim){
              $max_motif=$max_sim;
              $max_cap=$cap1;
              @max1=@max_pwm1;
              @max2=@max_pwm2;
          }
          $grade=$max_sim;
          print $cap2."\t".$cap1."\t".$grade."\n";
          #print "> ".$cap2."\t".$cap1."\t".$max_sim."\n";
          #showpwm(\@max_pwm1,0);
          #showpwm(\@max_pwm2,0);
          undef @max_pwm1;
          undef @max_pwm2;
        
     }
     #print $cap2."\t".$grade."\n";
     
     
     #print "> ".$cap2."\t".$max_cap."\t".$max_motif."\n";
     #showpwm(\@max1,0);
     #showpwm(\@max2,0);
     undef @max1;
     undef @max2;
}

