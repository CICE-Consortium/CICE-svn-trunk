#-----------------------------------------------------------------------
#  This script generates a makefile dependency list for sources
#  which occur in several directories.
#  It accepts a list of directories as arguments:
#    the first directory is the name of a compilation working directory
#    the remaining directories are locations for source code
#
#  author: Phil Jones, LANL
#  modified by Elizabeth Hunke, LANL
#-----------------------------------------------------------------------

use English;

#-----------------------------------------------------------------------
#  Create a working directory from the first argument
#-----------------------------------------------------------------------

$compdir = shift(@ARGV);
if (! -e $compdir) {
  print "Creating working directory $compdir for compilation\n";
  mkdir $compdir, 0775;
}

#-----------------------------------------------------------------------
#  Create Depends file for file dependencies and source file for list
#  of sources.
#-----------------------------------------------------------------------

open(DEPFILE, "+>Depends") or die "Unable to create Depends/n";

#-----------------------------------------------------------------------
#  Initialize list of object and source files
#-----------------------------------------------------------------------

$nfiles = 0;
$objects[0] = "OBJS = \\\n";
$sources[0] = "SRCS = \\\n";

#-----------------------------------------------------------------------
#  Read each directory
#-----------------------------------------------------------------------

foreach $dir (@ARGV) {

   if (! -d $dir) {
      print "$dir is not a directory;";
      next;
   }

   if (! opendir(DIR, $dir) ) {
      print "Cannot open $dir\n";
      next;
   }

   print "Reading directory $dir \n";

#-----------------------------------------------------------------------
#  Process each file in the directory
#-----------------------------------------------------------------------

   @files = readdir(DIR);
   fileloop: foreach $file (@files){

     ###
     ### Skip over parent, present directory entries
     ###

     if ($file eq ".")  {next fileloop;}
     if ($file eq "..") {next fileloop;}

     ###
     ### Only examine files with suffixes F,F90
     ### and extract root filename
     ###

     $srcfile = $file; # save current filename as source file
     if ($file=~ s/(\w).(F90|F)\Z/\1/ == 0) {next fileloop;}

     if ($srcfile eq "$file.F") 
       {$suffix    = "F";
        $newsuffix = "f";}
     if ($srcfile eq "$file.F90") 
       {$suffix    = "F90";
        $newsuffix = "f90";}
     
     ###
     ### Increment file count and extract root filename
     ###

     $nfiles++;  # increment total number of files

     print "Making dependencies for file $srcfile\n";

#-----------------------------------------------------------------------
#    Add file to lists of object and source files
#-----------------------------------------------------------------------

     $objects[$nfiles] = "\t$compdir/$file.o \\\n";
     $sources[$nfiles] = "\t$dir$srcfile \\\n";
     
#-----------------------------------------------------------------------
#    Read file to check for dependencies and create a dependency for
#    the source file.
#-----------------------------------------------------------------------

     open(INFILE, "<$dir/$srcfile") or die "Error opening $dir$srcfile\n";

     print DEPFILE "$compdir/$file.o: $compdir/$file.$newsuffix\n";
     print DEPFILE "$compdir/$file.$newsuffix: $srcfile\n";

     while (<INFILE>){

#-----------------------------------------------------------------------
#       Check to see if this line is an F90 use statement and
#       if it is, create a dependency entry
#-----------------------------------------------------------------------

        # check to see if this line is an F90 use statement
        if (m/(^\s*)(u|U)(s|S)(e|E)(\s+)(\w+)(\s|\S)*$/){

          # extract module name from a use line
          $modname = $ARG;
          $modname =~ s/(^\s*)(u|U)(s|S)(E|e)(\s+)(\w+)(\s|\S)*$/\6/;

          # create entry into dependency file
          print DEPFILE "$compdir/$file.o: $compdir/$modname.o\n";
        }
     }

     close(INFILE);
     
   }

   closedir(DIR);

}

#-----------------------------------------------------------------------
#  Create an Object file and write list of objects to this file
#-----------------------------------------------------------------------

###
### Remove continuation character from last line
###

$objects[$nfiles] =~ s/(o\s\\)/o/g;

open(OBJFILE, "+>Objects") or die "Unable to create Objects/n";
for ($i=0; $i<=$nfiles; $i++){
   print OBJFILE $objects[$i];
}
close(OBJFILE);

#-----------------------------------------------------------------------
#  Create a Sources file and write list of sources to this file
#-----------------------------------------------------------------------

###
### Remove continuation character from last line
###

$sources[$nfiles] =~ s/(\s\\)//g;

open(SRCFILE, "+>Sources") or die "Unable to create Sources/n";
for ($i=0; $i<=$nfiles; $i++){
   print SRCFILE $sources[$i];
}
close(SRCFILE);

#-----------------------------------------------------------------------
# All done
#-----------------------------------------------------------------------

