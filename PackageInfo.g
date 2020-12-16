#
# MotivesForBiArrangements: Orlik-Solomon bicomplexes for bi-arrangements
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "MotivesForBiArrangements",
Subtitle := "Orlik-Solomon bicomplexes for bi-arrangements",

Version := Maximum( [
                   "2017.05.01", ## Mohamed's version
                   ## this line prevents merge conflicts
                   "2017.05.01", ## Daniel's version
                   "2018.02.07"
                   ] ),

Date := ~.Version{[ 1 .. 10 ]},
Date := Concatenation( ~.Date{[ 9, 10 ]}, "/", ~.Date{[ 6, 7 ]}, "/", ~.Date{[ 1 .. 4 ]} ),

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Mohamed",
    LastName := "Barakat",
    WWWHome := "http://www.mathematik.uni-kl.de/~barakat/",
    Email := "mohamed.barakat@uni-siegen.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
  rec(
  	IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Clément",
    LastName := "Dupont",
    WWWHome := "http://www.math.univ-montp2.fr/~dupont/",
    Email := "clement.dupont@umontpellier.fr",
    PostalAddress := Concatenation(
            	"Institut Montpelliérain Alexander Grothendieck (IMAG)\n",
            	"Université de Montpellier\n",
            	"Place Eugène Bataillon\n",
            	"34090 Montpellier\n",
              	"France" ),
    Place := "Montpellier",
    Institution := "Université de Montpellier",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Daniel",
    LastName := "Juteau",
    WWWHome := "https://webusers.imj-prg.fr/~daniel.juteau/",
    Email := "daniel.juteau@imj-prg.fr",
    PostalAddress := Concatenation(
               "Institut de Mathématiques de Jussieu\n",
               "Paris Rive Gauche (IMJ-PRG)\n",
               "UP7D - Campus des Grands Moulins\n",
               "Boite Courrier 7012\n",
               "Bâtiment Sophie Germain\n",
               "8 Place Aurélie Nemours,\n",
               "75205 PARIS Cedex 13",
               "France" ),
    Place := "Paris",
    Institution := "Institut de Mathématiques de Jussieu - Paris Rive Gauche (IMJ-PRG)",
  ),
],

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/homalg-project/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
#SupportEmail   := "TODO",
PackageWWWHome  := "https://homalg-project.github.io/MotivesForBiArrangements/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "MotivesForBiArrangements",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Orlik-Solomon bicomplexes for bi-arrangements",
),

Dependencies := rec(
  GAP := ">= 4.8",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.5" ],
                   [ "alcove", ">= 2016-09-21" ],
                   [ "LinearAlgebraForCAP", ">= 2018.02.07" ],
                   [ "Bicomplexes", ">= 2017.05.11" ]
                   ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

Keywords := [ "Orlik Solomon", "bicomplexes" ],

));


