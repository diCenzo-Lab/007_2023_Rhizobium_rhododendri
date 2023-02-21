#!usr/bin/perl
use 5.010;

# Get file names
$motifFile = @ARGV[0];
$gffFiles = @ARGV[1];
$genomeFile = @ARGV[2];

# Modify gffFiles name in preparation for output files
$gffFilesOutput = $gffFiles;
$gffFilesOutput =~ s/_modified.gff//;
$gffFilesOutput =~ s/Data/Intermediate/;

# Extract replicons
open($in, '<', $gffFiles);
while(<$in>) {
  chomp;
  if(/sequence-region/) {
    @line = split(' ', $_);
    push(@replicons, @line[1]);
  }
}
close($in);

# Extract sequence of each replicon and then calculate GC skew
foreach $i (@replicons) {
  $test = 0;
  open($in, '<', $genomeFile);
  $intermediate = $gffFilesOutput . '_' . $i . '_' . 'sequence.fna';
  open($out, '>', $intermediate);
  while(<$in>) {
    if(/>/) {
      $_ =~ s/\./_/g;
      $_ =~ s/\ /_/g;
      $_ =~ s/\,//g;
      $_ =~ s/\(//g;
      $_ =~ s/\)//g;
      $_ =~ s/\-/_/g;
      if(/$i/) {
        $test = 1;
        print $out ($_);
      }
      else {
        $test = 0;
      }
    }
    else {
      if($test == 1) {
        print $out ($_);
      }
    }
  }
  close($in);
  close($out);
  open($in, '<', $intermediate);
  $output = $gffFilesOutput . '_' . $i . '_' . 'GCskew.txt';
  $output =~ s/Intermediate/Output/;
  open($out, '>', $output);
  while(<$in>) {
    if(/>/) {
      $G = 0;
      $C = 0;
      @data = ();
    }
    else {
      chomp;
      @line = split('', $_);
      foreach $j (@line) {
        push(@data, $j);
      }
    }
  }
  close($in);
  $size = scalar @data;
  for($n = 1; $n <= $size; $n = $n + 10000) {
    $G = 0;
    $C = 0;
    $min2 = $n - 9999;
    if($min2 < 1) {
      $min = $size + $min2;
    }
    else {
      $min = $min2;
    }
    if($min2 < 1) {
      for($a = 1; $a <= $n; $a++) {
        if(@data[$a] eq 'G') {
          $G++;
        }
        if(@data[$a] eq 'C') {
          $C++;
        }
      }
      for($a = $min; $a <= $size; $a++) {
        if(@data[$a] eq 'G') {
          $G++;
        }
        if(@data[$a] eq 'C') {
          $C++;
        }
      }
      $skew = ($G - $C) / ($G + $C);
      say $out ("$min\t$n\t$skew");
    }
    if($min2 >= 1) {
      for($a = $min; $a <= $n; $a++) {
        if(@data[$a] eq 'G') {
          $G++;
        }
        if(@data[$a] eq 'C') {
          $C++;
        }
      }
      $skew = ($G - $C) / ($G + $C);
      say $out ("$min\t$n\t$skew");
    }
  }
  $G = 0;
  $C = 0;
  for($a = $size - 9999; $a <= $size; $a++) {
    if(@data[$a] eq 'G') {
      $G++;
    }
    if(@data[$a] eq 'C') {
      $C++;
    }
  }
  $skew = ($G - $C) / ($G + $C);
  $min = $size - 9999;
  say $out ("$min\t$size\t$skew");
  close($out);
  $content = '';
  open($in, '<', $output);
  while(<$in>) {
    $content = $content . $_;
  }
  close($in);
  open($out, '>', $output);
  say $out ("Start_nt\tEnd_nt\tGC_skew");
  print $out ($content);
  close($out);

}
