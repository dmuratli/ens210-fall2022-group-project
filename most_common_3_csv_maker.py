iFile = open( "results1000.txt" , "r")
listo = []
for elem in iFile:
    listo.append( elem ) 
iFile.close()
most_L = []

for elem in range( len( listo ) ):
    if listo[elem] == "\n":

        num = listo[ elem - 5 ].split( " " )[-1][0:-1]
        
        most_L.append( num + "|" + listo[ elem - 15 ][1:-1] )
        most_L.append( num + "|" + listo[ elem - 14 ][1:-1] ) 
        most_L.append( num + "|" + listo[ elem - 13 ][1:-1] )
        


Cansu = open( "most_c_3.csv", "a")

for elem in range( len( most_L ) ):
    inTxt = most_L[elem].replace( " ", "|")
    inTxt = inTxt.replace("with|", "" )
    inTxt = inTxt.replace("|", ",")[0:-1]
    Cansu.write( inTxt )
    Cansu.write( "\n" )

Cansu.close()


