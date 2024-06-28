'''
Rizoan Azfar - CPS109
501237799
Program Description: DNA Transcriber and Translator.

This program reads a file which comprises of a strand of DNA (The entire file is treated as one strand, overlooking newlines and spaces). The program will determine whether to use the provided file or its complement as the coding strand, and transcribe it to RNA and then translate the RNA to an amino acid sequence, and store them in separate files. In the case of a missing start codon, the sequence will not be transcribed or translated, and the user will be informed through the terminal. In the case of missing stop codons, the sequence will be transcribed until the end of the provided sequence and alert the user. 
'''

import time

#dictionary storing the 20 amino acids and their corresponding codons
#global variable that can be used if this project is expanded in the future
aminoLib = {
    "Phe": ['UUU', 'UUC'],
    "Leu": ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    "Ile": ['AUU', 'AUC', 'AUA'],
    "Met": ['AUG'],
    "Val": ['GUU', 'GUC', 'GUA', 'GUG'],
    "Ser": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    "Pro": ['CCU', 'CCC', 'CCA', 'CCG'],
    "Thr": ['ACU', 'ACC', 'ACA', 'ACG'],
    "Ala": ['GCU', 'GCC', 'GCA', 'GCG'],
    "Tyr": ['UAU', 'UAC'],
    "His": ['CAU', 'CAC'],
    "Gln": ['CAA', 'CAG'],
    "Asn": ['AAU', 'AAC'],
    "Lys": ['AAA', 'AAG'],
    "Asp": ['GAU', 'GAC'],
    "Glu": ['GAA', 'GAG'],
    "Cys": ['UGU', 'UGC'],
    "Trp": ['UGG'],
    "Arg": ['CGU', 'CGC', 'CGA', 'CGG', 'AGG', 'AGA'],
    "Gly": ['GGU', 'GGC', 'GGA', 'GGG'],
    "stop": ['UAA', 'UAG', 'UGA']
}

def fileReader(fileName: str) -> list:
    """Attmepts to read specified file name and returns contents in a list.

        Parameters:
            fileName: Name of file desired to be read

        Returns:
            basePairs: List comprising of each base pair read from the file
    """

    try:
        file = open(fileName, 'r')
        basepairs = file.read()
        
    except Exception as e:
        raise e
    else:
        file.close()
        print("File:", fileName, "read successfully")
        return basepairs
        #file will be read and closed if there are no issues, and print to the user
    
def fileWriter (fileName: str, info: any, writeType:bool = 0 ) -> None:
    """Writes or appends provided information to be stored in a file.

        Parameters:
            fileName: Name of file desired to be written to
            info: The data that will be written
            writeType: Boolean, defaults 0 = overwrite, 1 = append.

        Returns:
            None, the function only writes to the specified file
    """

    if writeType != 0:
        openType = 'a'
    else:
        openType = 'w'
    #changing file write type if provided value is changed
    
    file = open(fileName, openType)


    if writeType == 1:
        file.write("Protein sequence generated at: " + time.strftime("%c, %Z") + "\n")
        #includes the time at which the sequence was generated

    file.write(str(info) + "\n")
    #writes info and adds extra line for formatting

    file.close()
    #closing file immediately after writing, saving resources
    print("File Writing to '" + fileName + "' successful")

def dnaComplement(strand: str) -> str:
    """Converts provided DNA strand into its complementary strand.
    """

    complement = ""
    #goes through each element in the strand and adds its complementary DNA base pair (A<->T, G<->C)
    for i in strand:
        if i == 'A':
            complement += 'T'
        elif i == 'T':
            complement += 'A'
        elif i == 'G':
            complement += 'C'
        elif i == 'C':
            complement += 'G'
        else:
            complement += ''
            #non DNA pairs are ignored in the complementary strand
            #if a DNA strand with errors already exists, it is unlikely for the parent strand to be corrected since it is not the strand being replicated. Existing errors can cause major problems in the synthesis of proteins. 

    return complement

def transcriptionStart(codStrand: str, compStrand: str) -> tuple:
    """Finds the optimal starting index for transcription.

        Parameters:
            codStrand: Initial DNA strand
            compStrand: Complementary DNA strand

        Returns:
            Size 2 tuple where the first index is an int of the index for transcription, and the second index is a bool where False/0 will signify initial and True/1 will signify the complementary strand.
    """

    if 'ATG' not in codStrand and 'ATG' not in compStrand:
        print("There is no adequate start codon, DNA transcription cannot occur.")
        return None
    #DNA transcription requires the presence of ATG to exist to create the AUG from the coding strand/non template strannd. If ATG is not detected in either provided or complementary strand, RNA translation will not occur.

    elif 'ATG' in codStrand:
        codStart = codStrand.index('ATG')
        #store the earliest index at a possible starting point for transcription
        Start = codStart
    
        if 'ATG' in compStrand:
            #only occurs if both strands contain ATG

            compStart = compStrand.index('ATG')
            #stores the index of the starting point in the other strand


            #conditionals to check the earlier instance between the two strands to use as the non-template strand. Returns the index of the earlier strand as well as a boolean to represent which strand is to be used
            if codStart < compStart:
                return (Start, 0)
            else:
                Start = compStart
                return (Start, 1)
                
        else:
            return(Start,0)
        #returns index of the initial strand if ATG doesnt exist in both

    elif 'ATG' in compStrand:
        Start = compStrand.index('ATG')
        return (Start, 1)
    #returns index of complementary strand if ATG only exists in it.


def rnaStrand(codStrand: str, compStrand: str) -> str:
    """Transcribes DNA into RNA.

        Parameters:
            codStrand: Initial DNA strand
            compStrand: Complementary DNA strand

        Returns:
            rna: Resulting RNA strand after translating the optimal DNA strand
    """

    tInfo = transcriptionStart(codStrand, compStrand)
    #obtains the starting index and which strand to use

    if tInfo == None:
        raise Exception("Translation cannot occur")
    #Since no start codon was found using the transcriptionStart function, there will be no place for transcription to start

    rna = ""
    
    index = tInfo[0]
    #obtaining starting index from tInfo

    if tInfo[1] == 0:
        strand = codStrand
    else:
        strand = compStrand
    #conditional to determine which strand to use based on second element in tInfo

    while index < len(strand)-3:
    #using a while loop since transcription will continue until a stop codon is encountered, 
        
        for i in range(3):
            if strand[index+i] == 'T':
                rna += 'U'
            else:
                rna += strand[index+i]
            #The chosen strand is treated as the coding strand or the non template strand, meaning the resulting RNA sequence will be the same excluding Thymine, which doesnt exist in RNA and appears as Uracil
        index += 3
        #using a for loop to look at 3 codons at a time, as per what occurs in DNA transcription and translation. This means that the resulting RNA strand will always have a length divisble by 3

        endCodon = rna[len(rna)-3:]
        #looking at the current last 3 codons

        if endCodon == 'UGA' or endCodon =='UAA' or endCodon =='UAG':
            break
        #transcription will end if a stop codon is reached

        if index >= len(strand)-3:
            print("No stop codon found, improper protein may be synthesized incorrectly")
    
    print("DNA transcribed successfully")
    return rna
    
def translation(rna: str) -> list:
    """Translates RNA sequence into a sequence of amino acids.

        Parameters:
            rna: RNA strand to be translated

        Returns:
            aminoacids: list containing the translated sequence of amino acids.
    """
    aminoacids = []
    
    for i in range(0, len(rna), 3):
        #looking at 3 codons at once

        #this loop iterates through the global dictionary at the top of the program. It checks for which amino acid the current codon will code for and appends it to the list that will be returned
        for key in aminoLib:
            if rna[i:i+3] in aminoLib[key]:
                aminoacids.append(key)
    print("RNA translated successfully")
    return aminoacids


if __name__ == "__main__":
    startTime = time.time()
    #will be used to calculate runtime

    codingStrand = fileReader("a1_inputDNA.txt")
    complementaryStrand = dnaComplement(codingStrand)
    #reading the specified file which contains a strand of DNA that can be modified. The complementary DNA strand is found immediately after reading it. 
    

    rnaFileName = "RNAOutput.txt"
    proteinFileName = "ProteinOutput.txt"
    #file names of the rna and protein files, can be modified if needed
    
    rna = rnaStrand(codingStrand, complementaryStrand)
    #obtain the RNA strand
    fileWriter(rnaFileName, rna)
    #rnaRead = fileReader(rnaFileName)
    #^^^RNA output file can be read. would be useful in a multi-file program

    proteins = translation(rna)
    fileWriter(proteinFileName, proteins, 1)
    #obtains the amino acid sequence and writes it to the protein output file


    endTime = time.time()
    timeElapsed = endTime-startTime
    #calculate the elapsed time taken to run the program

    print("\nDNA transcribed and translated in", "{:.4f}".format(timeElapsed), "seconds.\nPlease check", rnaFileName, "and", proteinFileName, "for your respective RNA and Protein outputs")
    #inform the user of the elapsed time, and file names for the results of transcription and translation