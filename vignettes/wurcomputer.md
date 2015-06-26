## Virtual machine

You can use virtual machine to access the chip databases. Ask admins for details.


## Database computer

To access HITChip database from WUR database computer, do the following:

1. Login to the new database machine. It will ask two passwords during login. You can get these from the [admins](contact)

1. Create a new output folder for profiling files. To do this, open the folder icon on the left panel, navigate to the "Users" directory and Create New Folder with right-mouse click.

1. Open terminal window (the black screen icon on the left panel).

1. Start R by giving the standard command 'R' in the terminal window.


### Extracting HITChip data from the MySQL database ('profiling')

Run the following commands in R (on a computer or virtual machine that
has access to the database). 

You can request the username, password, database name and port from
the [admins](contact). You also need SSH port forwarding, ask for more
instructions.

The program will guide you through preprocessing by asking to specify
the projects, samples, and preprocessing parameters. For further
instructions, contact the [admins](contact).



```r
# Load the package
library(HITChipDB)

# Run profiling script over FTP
# NOTE: need to ssh port forwarding for localhost first - detailed instructions from the admins!
# NOTE: request the correct dbpwd and port from the admins (see above)
params <- run.profiling.script(dbuser = "root", dbpwd = "array", dbname = "Phyloarray", host = '127.0.0.1', port = ask.correct.port.number.from.admins) 
```

This prompts you through profiling options described in [full instructions](Protocol.md).


### HITChip

To extract HITChip data, see [these instructions](Hitchip).


### MITChip

To extract MITChip data, give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "mit", dbpwd = "passu", dbname = "phyloarray_mit")
```

This prompts you through profiling options described in [full
instructions](protocol).


### PITChip (older chip)

To extract PITChip data (older chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "pit", dbpwd = "passu", dbname = "phyloarray_pit")
```

This prompts you through profiling options described in [full
instructions](protocol).



### PITChip (new chip)

To extract PITChip data (new chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "pit", dbpwd = "passu", dbname = "pitchipdb")
```

This prompts you through profiling options described in [full
instructions](protocol).


### ChickChip (older chip)

To extract ChickChip data (older chip), give the command in R:


```r
library(HITChipDB) 
params <- run.profiling.script(dbuser = "mit", dbpwd = "passu", dbname = "chickchipdb")
```

This prompts you through profiling options described in [full instructions](Protocol).




