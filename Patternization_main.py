#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:41:22 2024

@author: alexanderpfaff
"""

import json
from itertools import permutations, combinations
from collections import Counter
from random import choice 
from math import factorial
from copy import deepcopy



class Patternize:
    """
    Patternization: version NPEGL_0.9 
    
    WELCOME to >> P A T T E R N I Z A T I O N <<  and
    CONGRATULATIONS on your acquisition of this cutting-edge devise! 
    
    Patternization is a method for classifying, analysing and visualising
    syntactic / word-order variation on the basis of annotated text corpora.
    The current version immediately builds upon the NP database NPEGL.
    @see: Pfaff (2019)
    @see: Pfaff (2024)
    @see: Pfaff & Bouma (2024)
    
    """
    
    def __init__(self, DB_from_file='oice', _secondaryInput=None):
        """
        
        Initializes a Patternize::Database obkect comprising 
            -- a database (list) of NPs (type: Pattern),  and
            -- methods/functionalities to measure, process and classify 
                       syntactic / word order variation found in the database.
        This constructor provides two ways of initializing: 
            -- via primary input -- database is read in from file,
            -- via secondary input; this option is only used internally in 
                order to activate Patternize objects acting as attributes of 
                the parent object.

        Parameters
        ----------
        DB_from_file : 
            TYPE:        json file,  
            DESCRIPTION: the database to work on; 
                         the options are: 
                             - 'oice', 
                             - 'osax', 
                             - 'oeng'. 
                         The default is 'oice'.
                         
        _secondaryInput :   Internal use only!
            TYPE:           None or database.
            DESCRIPTION:    The default is None.

        Returns
        -------
        None.

        """
        if type(_secondaryInput) == list:
            self._database = _secondaryInput[1:]
            self._ID = _secondaryInput[0]
        else:
            with open(DB_from_file + "_db.json","r") as f:
                inputDatabase = json.load(f)
            self._database = [Pattern(np)
                              for np in inputDatabase]
            self._ID = inputDatabase[0]["Language"]
        
        self._maxPatternSize = None
        self._update = []
        self._rnd_update = 0
        self._rnd_database = None
        self._prenominalDomain = None
        self._postnominalDomain = None


###############################
        
    @property
    def ID(self) -> str:
        """
        Returns:    str; 
        -------     language of the current DB.
        
        """
        return self._ID

    @property
    def size(self):
        """
        Returns:    int; 
        -------     number of NPs in the current DB.
        
        """
        return len(self.database)
        
    @property
    def maxPatternSize(self):
        """
        Returns:    int; 
                    maximal pattern size = longest pattern in current DB;
        -------     measured in number of category labels.
        
        """
        return len(max(self.patternize(), key=len))
    
    @property
    def update(self):
        """
        Returns:     list; 
                     list of features according to which DB has been updated;
        -------
        DESCRIPTION: (X , bool): NPs containing X have been
                                   - retained (True), else removed (if absent);
                                   - removed (False), else retained (if absent).
                                     
        """
        return tuple(self._update)

    @property
    def rnd_update(self):
        """
        Returns:    int; 
                    number of random databases yet generated
        -------         (= number of calls >>> self.randomize())
        
        """
        return self._rnd_update

    @property
    def database(self):
        """
        THIS database containing all current NPs in >>> self.ID 
        
        Returns
        -------
        TYPE:       list; 
                    TYPE of elements in list: Pattern.
        
        """
        return self._database

    @property
    def rnd_database(self):
        """
        This recursive attribute is itself a Patternize database with all 
        functionalities. It contains a random selection of a given size of the
        mother DB. 
        It is initially empty, and can be activated by the method 
        >>> self.randomize(size)

        Returns
        -------
        TYPE:       Patternize
        
        """
        if self._rnd_database == None:
            self._rnd_database = Patternize("empty")
        return self._rnd_database

    
    
    @property 
    def prenominalDomain(self):
        """
        This recursive attribute is itself a Patternize database with all 
        functionalities. It contains the prenominal categories to the exclusion
        of the noun itself. 
        It is initially empty, and can be activated by the method 
        >>> self.partitionize(partition="pre")

        Returns
        -------
        TYPE:       Patternize
        
        """
        if self._prenominalDomain == None:
            self._prenominalDomain = Patternize("empty")
        return self._prenominalDomain

    @property 
    def postnominalDomain(self):
        """
        This recursive attribute is itself a Patternize database with all 
        functionalities. It contains the postnominal categories to the exclusion
        of the noun itself. 
        It is initially empty, and can be activated by the method 
        >>> self.partitionize(partition="post")

        Returns
        -------
        TYPE:       Patternize
        
        """
        if self._postnominalDomain == None:
            self._postnominalDomain = Patternize("empty")
        return self._postnominalDomain

    
    
###############################

        
    def filter(self, cat, present=True):
        """
        Allows to prune the NP entries in the current working database. 

        Parameters
        ----------
        cat : TYPE:         str (= category label),    
            DESCRIPTION:    the target category;
        present : TYPE:     bool; 
            DESCRIPTION:    True: NPs with target category will be removed,
                            False: NPs witout target category will be removed;
                             The default is True.

        Returns    
        -------
        None.

        """
        eliminate = []
        for np in self.database:
            if np.hasCat(cat) == present:
                eliminate.append(np)
        _ = [self.database.remove(np)
             for np in eliminate]
        self._update.append((cat, present))
        
        
    def ndb(self, properN=False, onlyProper=False, includeConjuncts=False):
        """
        'ndb' = n_ominal d_ata b_ase. 
        Creates a working DB by imposing two restrictions on THIS DB:
            -- every entry contains a noun ("N.C", "N.P" or "N"),
                --> default noun == "N.C" (= common noun)
            -- coordinated structures are excluded.
            @see: self.filter({noun}, False) & self.filter("&", True) 
        It is possible to include the individual cojuncts before 
        the coordination structure is discarded. 

        Parameters
        ----------
        properN :   
            TYPE:           bool, optional.
            DESCRIPTION:    Include proper nouns. The default is False.
            
        onlyProper : 
            TYPE:           bool, optional.
            DESCRIPTION:    Proper nouns only. The default is False.
            
        includeConjuncts : 
            TYPE:           bool, optional.
            DESCRIPTION.    Include the individual conjuncts into the DB 
                            BEFORE coordination is eliminated. 
                            The default is False.
        Returns
        -------
        None.

        """
        if includeConjuncts:
            self.partitionize(partition="pre", catCondition="&")
            self.partitionize(partition="post", catCondition="&")            
        noun = self._catCondition(catCondition="N.C", properN=properN, 
                                 onlyProper=onlyProper)
        if not (((noun, False) in self.update) and (("&", True) in self.update)):
            self.filter(noun, False)
            self.filter("&", True)
            self._update.append(("includeConjuncts", includeConjuncts))


################################

    def getNP(self, index): 
        """
        

        Parameters
        ----------
        index : TYPE:       integer
            DESCRIPTION:    index position in this database

        Returns:            the NP @{index}
        -------
        """
        return self.database[index]


    def patternize(self, level=2):
        """
        Pattern inventory of the current DB at level={level}.
        A 'pattern' is an NP type (not token),
        defined by the linear order of the categories involved.
        De facto a getter method. 
        
        Parameters
        ----------
        level : TYPE:        int; 
                DESCRIPTION: the (sub-)categorization level: 0-3, 7.
                             The default is 2. 

        Returns
        -------
                TYPE:        Counter; 
                DESCRIPTION: pattern : occurrences in current DB 
                             (=key)     (= value)

        """
        out = [patt.getPatt(level)
               for patt in self.database]
        return Counter(out)
    

    def categorize(self, level=2):
        """
        Category inventory of the current DB at level={level}.
        De facto a getter method. 

        Parameters
        ----------
        level : TYPE:        int;
                DESCRIPTION: the (sub-)categorization level: 0-3, 7;
                             The default is 2. 

        Returns
        -------
                TYPE:        Counter; 
                DESCRIPTION: category : occurrences in current DB 
                             (=key)     (= value)

        """
        out = [cat
               for np in self.database
               for cat in np.getPatt(level)]    ## instead of getCat !!!
        return Counter(out)                     ## --> multiple occurrences

    def getAllCats(self):
        """
        Category inventory of the current DB at all levels.
        De facto a getter method. 
        
        Returns
        -------
                TYPE:        Counter; 
                DESCRIPTION: category : occurrences in current DB 
                             (=key)     (= value)

        """
        out = [cat
               for np in self.database
               for cat in np.allCats]
        return Counter(out)
        



    def lemmatize(self):
        """
        Returns
        -------
        TYPE:             Counter;
            DESCRIPTION:  lemma found in current DB : occurrences in current DB
                          (= key)                     (= value)

        """
        out = [lemma.lower()
               for np in self.database
               for lemma in np.lemmata]
        return Counter(out)


    def findLemma(self, cat):
        """
        Finds all lemmata in THIS database instantiating {cat}.
        
        Parameters
        ----------
        cat : TYPE:       str (category label);

        Returns:          
        -------           
        TYPE:             Counter; 
            DESCRIPTION:  lemma instantiating {cat} : occurrences in current DB
                          (= key)                     (= value)

        """
        out = []
        for np in self.database:
            if cat in np.allCats:
                for position in range(np.length):
                    if cat in np.macroPattern[position]:
                        out.append(np.lemmata[position])
        return Counter(out)                    


    def cat_in_patt(self, cat, level=2):
        """
        returns the patterns in which a given category {cat} occurs.
        
        Parameters
        ----------
        cat : TYPE:         str (category label);
            DESCRIPTION.
        level : TYPE:       int;  
            DESCRIPTION:    level of sub-categorization.
                            The default is 2.

        Returns
        -------
        TYPE:               Counter;  
            DESCRIPTION:    pattern containing {cat} : occurrences in curent DB
                            (= key)                    (= value)

        """
        return Counter([np.getPatt(level)
                        for np in self.database
                        if cat in np.getPatt(level)])

        
    @property
    def tokenize(self):
        """
        Returns
        -------
        TYPE:             int;
            DESCRIPTION:  number of tokens in the current DB. 

        """
        return sum([np.length
                    for np in self.database])
    
    @property
    def lemmata_per_NPs(self):
        """
        Nomen est omen!

        Returns
        -------
        TYPE:             float; 
            DESCRIPTION:  ratio lemmata per NPs in the current DB.

        """
        return round(100 * (len(self.lemmatize()) / self.size), 2)
        
    @property
    def lemmata_per_tokens(self):
        """
        Nomen est omen!

        Returns
        -------
        TYPE:             float; 
            DESCRIPTION:  ratio lemmata per tokens in the current DB.

        """
        return round(100 * (len(self.lemmatize()) / self.tokenize), 2)
        
    def patterns_per_NPs(self, level=2):
        """
        Nomen est omen!

        Parameters
        ----------
        level : TYPE:     int;
            DESCRIPTION:  level of sub-categorization.
                          The default is 2.
             
        Returns
        -------
        TYPE:             float; 
            DESCRIPTION:  ratio patterns per NPs in the current DB.

           @see patternDiversity! 
        """
        return round(100 * (len(self.patternize(level)) / self.size), 2)
    
    def categories_per_patterns(self, level=2):
        """
        Nomen est omen!

        Parameters
        ----------
        level : TYPE:     int;
            DESCRIPTION:  level of sub-categorization.
                          The default is 2.

        Returns
        -------
        TYPE:             float;
            DESCRIPTION:  ratio categories per NPs in the current DB.

        """
        return round(100 * (len(self.categorize(level)) / 
                            len(self.patternize(level))), 2)



    def randomize(self, size=1000):
        """
        Creates a randomized database of a given size {size} == 
        subset of the mother database. 

        Parameters
        ----------
        size : TYPE:       int;
            DESCRIPTION:   size of the random sub-database.  
                           The default is 1000.

        Returns
        -------
        None.

        """
        out = set()
        while len(out) < size:
            out.add(choice(self.database))
        out = list(out)
        out.insert(0, self.ID + ", @rnd_database")
                
        self._rnd_update +=1
        self._rnd_database = Patternize(DB_from_file=None, _secondaryInput=out)


    def patternDiversity(self, level=2, runs=500, size=1000): 
        """
        Pattern Diversity = ratio patterns at {level} per NPs relative to a 
                             standardized common denominator. 
                             The ratio is determined as the means of {runs} 
                             iterations of patterns per NPs over randomized 
                             sub-DBs of size {size}.

        Parameters
        ----------
        level : TYPE:     int;
            DESCRIPTION:  level of sub-categorization.
                          The default is 2.
        runs : TYPE:      int;  
            DESCRIPTION:  number of iterations over random sub-DBs.
                          The default is 500.
        size : TYPE:      int; 
                          size of the random sub-DB.
            DESCRIPTION.  The default is 1000.

        Returns
        -------
        TYPE:             float; 
            DESCRIPTION:  pattern diversity (ratio: patterns per {size} NPs)

        """
        patterns_from_subDBs = 0
        for run in range(runs):
            self.randomize(size=size)
            patterns_from_subDBs +=len(self.rnd_database.patternize(level=level))
        avg_patts = patterns_from_subDBs / runs
        return round((100 * (avg_patts / size)), 2)
            
    
    def clear_subDBs(self):
        """
        Sets all activated sub-databases to "empty".

        Returns
        -------
        None.

        """
        self._rnd_database = Patternize("empty")
        self._prenominalDomain = Patternize("empty")
        self._postnominalDomain = Patternize("empty")
    
    

###############################################
#########    COMBINATORIX


    def _catCondition(self, catCondition="N.C", properN=False, 
                     onlyProper=False):
        """
        auxiliary method        

        Parameters
        ----------
        catCondition : TYPE: str or a neg-type (None/False);
            DESCRIPTION:     category label.
                             The default is "N.C".
        properN : TYPE:      bool;  
            DESCRIPTION:     include proper names. 
                             The default is False.
        onlyProper : TYPE:   bool; 
            DESCRIPTION:     proper names only.
                             The default is False.

        Returns
        -------
        catCondition : TYPE: str or a neg-type (None/False)
            DESCRIPTION:     category label to be used in filter methods. 

        """
        if catCondition:
            if properN:
                if onlyProper:
                    catCondition = "N.P"
                else:
                    catCondition = "N"
        return catCondition


    def combinatorialFlexibility(self, cats, length=3, sm_pattern=None, 
                                 align=None, pattern_threshold=1, 
                                 group_threshold=2, count=bool, 
                                 addAdjectives=False, catCondition="N.C",
                                 properN=False, onlyProper=False):
        """
        This method probes into the range of combinatorial flexibility of the
        current DB with respect to the pattern permutations of a given size
        created on the basis of a given selection of category labels. 
        
        Procedure: with n = len(cats) and k = {length}, all n_C_k combinations
                   c in {cats} are produced.  If not catCondition or 
                   catCondition in c, all k! permutations p of c are generated 
                   and collected as permutation groups.  For each p in c, the
                   current DB is browsed and every match according to the
                   criterion sm_pattern is counted; if number of matches p >=
                   pattern_threshold, and added number of matches p in c >= 
                   group_threshold, the permutation group c is considered 
                   attested and stored as a CombFlex object.
                   Attestation values can be
                   - categorical (yes/no) <==> count=bool, or 
                   - numerical (number of attestations) <==> count=int. 

        Parameters
        ----------
        cats : TYPE:        container (list, tuple) of category labels;
            DESCRIPTION:    sample space on basis of which combinations of 
                            length {length} are created.
        length : TYPE:      int;  
            DESCRIPTION.    length of combinations/permutations.
                            The default is 3.
        sm_pattern : TYPE:  container (list, tuple) of category labels;
            DESCRIPTION:    search/match pattern template,
                            possible specification: >>> self.precise,
                                                    >>> self.rigid,
                                                    >>> self.flexi.
                            The default is None ==> rigid.
        align : TYPE:       str, optional; 
            DESCRIPTION:    allows specification for alignment; 
                            possible values: "left", "right".
                            The default is None.
        pattern_threshold : 
            TYPE:           int;  
            DESCRIPTION:    minimal number of attestations in order for a given 
                            permutation to be counted.
                            The default is 1.
        group_threshold : 
            TYPE:           int;
            DESCRIPTION:    minimal number of attestations of permutations 
                            within a given permutation group in order for the
                            combination/permutation group to be counted.
                            The default is 2.
        count : TYPE:       <int> or <bool>; 
            DESCRIPTION:    specifies the mode of counting: 
                            -- number of attestations (int),
                            -- True/False if (not) attested 
                            in accordance with the threshold settings.
                            The default is bool.        
        addAdjectives : 
            TYPE:           bool, optional; 
            DESCRIPTION.    if True, adjective cat labels will be added again;
                            ensures combinations --> permutations involving 
                            two adjectives are taken into account.
                            The default is False.
        catCondition : 
            TYPE:           str (category label);  
            DESCRIPTION:    specifies whether a given category must be present
                            in all combinations (--> permutations).
                            The default is "N.C".
        properN : TYPE:     bool, optional; 
            DESCRIPTION:    include proper names.
                            The default is False.
        onlyProper : TYPE:  bool, optional; 
            DESCRIPTION:    proper names only.
                            The default is False.

        Returns
        -------
        out : TYPE:         list [ {CombFlex} ]; 
            DESCRIPTION:    permutation groups satisfying the matching and 
                            threshold conditions are construed as CombFlex 
                            objects and appended to {out}.

        """
        catCondition = self._catCondition(catCondition=catCondition, 
                                          properN=properN, 
                                          onlyProper=onlyProper)        
        if sm_pattern == None:
            sm_pattern = self.rigid             
        out = []
        possiblePermutations = self._possiblePermutations(cats, length=length, 
                                                          catCondition=
                                                          catCondition,
                                                          addAdjectives=
                                                          addAdjectives)
        for combination in possiblePermutations:
            permutationGroup =  possiblePermutations[combination]
            groupAttestation = dict()
            groupCount = 0
            for i in range(factorial(length)):
                pattern = next(permutationGroup)
                attested = self.collectCustomMatch(pattern, 
                                                   sm_pattern=sm_pattern, 
                                                   align=align)
                if len(attested) < pattern_threshold:
                    attested = []
                groupAttestation[pattern] = dict()
                groupAttestation[pattern]["Count"] = count(len(attested))
                groupAttestation[pattern]["Attestations"] = attested
                groupCount += len(attested)
            if groupCount >= group_threshold:                
                combFlexOut = CombFlex(combination=combination, length=length, 
                                       sm_pattern=str(sm_pattern)[:-25] + ">", 
                                       align=str(align), 
                                       pattern_threshold=pattern_threshold, 
                                       group_threshold=group_threshold, 
                                       count=count, groupCount=groupCount,
                                       permutations=groupAttestation,
                                       catCondition=catCondition)
                out.append(combFlexOut)
        return out
    
        
 
    def customize(self, np, pattern, sm_pattern=None, align=None): 
        """
        Checks whether the NP matches a given pattern with additional 
        specifications provided by the functional parameter sm_pattern, which 
        can be set to self.precise, self.rigid, self.flexi (default: rigid);
        @see: precise, rigid, flexi. 
        In addition, it allows to impose an aligment condition with the 
        possible values "left": first category of m_pattern and s_pattern must
        be identical; "right":  last category of m_pattern and s_pattern must
        be identical; default: no alignment.  

        Parameters
        ----------
        np : TYPE:          Pattern;
            DESCRIPTION:    an NP. 
        pattern : TYPE:     list or tuple
            DESCRIPTION.
        sm_pattern : TYPE:  container (list, tuple) of category labels;
            DESCRIPTION:    search/match pattern template,
                            possible specification: >>> self.precise,
                                                    >>> self.rigid,
                                                    >>> self.flexi.
                            The default is None ==> rigid.
        align : TYPE:       str, optional; 
            DESCRIPTION:    allows specification for alignment; 
                            possible values: "left", "right".
                            The default is None.

        Returns
        -------
        TYPE:               bool;
            DESCRIPTION:    True if alignment and matching conditions are 
                            satisfied.

        """
        if sm_pattern == None:
            sm_pattern = self.rigid
        m_pattern = np.macroPattern
        s_pattern = list(pattern)
        # s_pattern = [cat for cat in pattern]
        if align == "left":
            alignmentCheck = s_pattern[0] in m_pattern[0]
        elif align == "right":
            alignmentCheck = s_pattern[-1] in m_pattern[-1]
        else: 
            alignmentCheck = True
        return alignmentCheck and sm_pattern(m_pattern, s_pattern)
        


    def precise(self, m_pattern, s_pattern):
        """
        A functional parameter ("pseudo-lambda") specifying a precise search 
        pattern (s_pattern) that has to match the given NP pattern (m_pattern). 
        
        Matching conditions:
        (Cat_1, Cat_2 .. Cat_n)  <==> (Cat_1, Cat_2 .. Cat_n)         
        i.e. an exact match.
        
        To be passed in to >>> self.combinatorialFlexibility as a 
        specification of the parameter 'sm_pattern';
        @see: flexi
        @see: rigid

        Parameters
        ----------
        m_pattern :        the matched pattern (NP)

        s_pattern :        the search pattern

        Returns
        -------
        TYPE:              bool;
            DESCRIPTION:   True if the matching condition is satisfied.
        """
        if len(m_pattern) != len(s_pattern):
            return False
        return all([s_pattern[index] in m_pattern[index]
                    for index in range(len(s_pattern))])


    def rigid(self, m_pattern, s_pattern):   
        """
        A functional parameter ("pseudo-lambda") specifying a rigid search 
        pattern (s_pattern) that has to match the given NP pattern (m_pattern). 
        
        Matching conditions:
        (Cat_1, Cat_2 .. Cat_n)  <==> (... Cat_1, Cat_2 .. Cat_n ...)         
        i.e. the sequence of s_pattern must be contained in m_pattern.
        
        To be passed in to >>> self.combinatorialFlexibility as a 
        specification of the parameter 'sm_pattern';
        @see: flexi
        @see: precise

        Parameters
        ----------
        m_pattern :        the matched pattern (NP)

        s_pattern :        the search pattern

        Returns
        -------
        TYPE:              bool;
            DESCRIPTION:   True if the matching condition is satisfied.
        """
        if len(s_pattern) > len(m_pattern):
            return False
        firstCatIndex = -1
        for index in range(len(m_pattern)):
            if s_pattern[0] in m_pattern[index]:
                firstCatIndex = index
        return self.precise(m_pattern[firstCatIndex:], s_pattern)
            

    def flexi(self, m_pattern, s_pattern):   
        """
        A functional parameter ("pseudo-lambda") specifying a flexible search 
        pattern (s_pattern) that has to match the given NP pattern (m_pattern). 
        
        Matching conditions:
        (Cat_1, Cat_2 .. Cat_n)  <==> (... Cat_1, ... Cat_2 ... Cat_n ...)         
        i.e. the relative sequence of categories in s_pattern must be found in
        m_pattern (regardless of intervening material).
        
        To be passed in to >>> self.combinatorialFlexibility as a 
        specification of the parameter 'sm_pattern';
        @see: precise
        @see: rigid

        Parameters
        ----------
        m_pattern :        the matched pattern (NP)

        s_pattern :        the search pattern

        Returns
        -------
        TYPE:              bool;
            DESCRIPTION:   True if the matching condition is satisfied.
        """
        if len(s_pattern) > len(m_pattern):
            return False
        matches = len(s_pattern)
        s_pattern.append("")
        match = False
        currentCat = s_pattern.pop(0)
        currentIndex = 0
        while currentIndex < len(m_pattern) and not match:
            if currentCat in m_pattern[currentIndex]:
                matches -=1
                if matches == 0:
                    match = True
                currentCat = s_pattern.pop(0)
            currentIndex +=1
        return match
                
            
        

    def collectCustomMatch(self, pattern, sm_pattern=None, align=None):
        """
        Collects NPs that satisfy the matching and alignment conditions.

        Parameters
        ----------
        pattern : TYPE:     container (tuple/list) of str; 
            DESCRIPTION:    sequence of category labels
        sm_pattern : TYPE:  functional parameter specifying matching conditions
            DESCRIPTION:    possible values: self.precise, self.rigid, 
                                             self.flexi
                            The default is None ==> rigid.
        align : TYPE:       str, optional;
            DESCRIPTION:    specifies an alignment condition: "left" or "right".
                            The default is None (= no alignment).

        Returns
        -------
        list                of Pattern objects;
            DESCRIPTION:    contains the NPs satisfying the matching and 
                            alignment conditions.

        """
        return [np
                for np in self.database
                if self.customize(np, pattern, sm_pattern=sm_pattern, align=align)]





    def combis_mirror(self, cats, sm):
        """ TO DO"""
        pass

    def combis_flanked(self, cats, sm):
        """ TO DO"""
        pass



    def _permutePatterns(self, pattern):
        """
        aux_@_combinatorialFlexibility
        creates all len(pattern)! permutations of {pattern}

        Parameters
        ----------
        pattern : TYPE:   tuple/list
            DESCRIPTION:  sequence of category labels

        Raises
        ------
        ValueError
            DESCRIPTION:  if len(pattern) < 2

        Returns
        -------
        TYPE:             iterator 
            DESCRIPTION:  contains all permutations of {pattern}.

        """
        if len(pattern) < 2:
            raise ValueError("Inappropriate size values!")
        return permutations(pattern)


    def _catCombinations(self, cats, length=3, catCondition="N.C"):
        """
        aux_@_combinatorialFlexibility
        creates all n_C_k combinations c over {cats}, with n = len(cats),
                                                                k = {length};
        if not catCondition or catCondition in c, c is added to output list.

        Parameters
        ----------
        cats : TYPE:      str; 
            DESCRIPTION:  category label
        length : TYPE:    int;
            DESCRIPTION:  length of combination (subset) in {cats}.
                          The default is 3.
        catCondition : 
            TYPE:         None or str, optional;
            DESCRIPTION.  specifies a category that must be contained in every
                          combination; normally a nominal category.
                          The default is "N.C".

        Raises
        ------
        ValueError
            DESCRIPTION:  if len(cats) < 2 || length < 2 || len(cats) < length.


        Returns
        -------
        TYPE:             tuple
            DESCRIPTION:  of tuples = combinations of length {length},
                          potentially satifying {catCondition}.

        """
        if len(cats) < 2 or length < 2 or len(cats) < length:
            raise ValueError("Inappropriate size values!")
        out = []
        combi_nations = combinations(cats, r = length)
        try:
            while True:
                combi = next(combi_nations)
                if not catCondition or (catCondition in combi):
                    out.append(combi)
        except StopIteration:
            return tuple(out) 


    def _possiblePermutations(self, cats, length=3, catCondition="N.C",
                              addAdjectives=False):
        """
        aux_@_combinatorialFlexibility

        Parameters
        ----------
        cats : TYPE:      str; 
            DESCRIPTION:  category label
        length : TYPE:    int;
            DESCRIPTION:  length of combination (subset) in {cats}.
                          The default is 3.
        catCondition : 
            TYPE:         None or str, optional;
            DESCRIPTION.  specifies a category that must be contained in every
                          combination; normally a nominal category.
                          The default is "N.C".
        addAdjectives : 
            TYPE:         bool, optional; 
            DESCRIPTION.  if True, adjective cat labels will be added again;
                          ensures combinations --> permutations involving 
                          two adjectives are taken into account.
                          The default is False.

        Returns
        -------
        out : TYPE:       dict: tuple of cat labels  :  iterator
            DESCRIPTION:        combination(s)       :  permutations
                                (= keys)             :  (= values). 

        """
        out = dict()
        cats = list(cats)
        if addAdjectives: 
            cats += ["Md.Aj", "Md.Aj.Fn", "Md.Aj.Lx"]
        for combi in self._catCombinations(cats, length=length, 
                                           catCondition=catCondition):
            out[combi] = self._permutePatterns(combi)
        return out
        

###############################################
###
###     partition the NP

    
    def partitionize(self, partition="pre", catCondition="N.C", properN=False, 
                     onlyProper=False):
        """
        Initializes NP partitions accessible via 
            --> self.prenominalDomain,
            --> self.postnominalDomain.
        Partitions are themselves Patternize objects

        Parameters
        ----------
        partition : TYPE: str, optional;
                    possible values: "pre", "post".
            DESCRIPTION. The default is "pre".
        catCondition : str, optional
            DESCRIPTION. The default is "N.C".
        properN : str, optional
            DESCRIPTION. The default is False.
        onlyProper : str, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        catCondition = self._catCondition(catCondition=catCondition, 
                                         properN=properN, 
                                         onlyProper=onlyProper)
        subDB = deepcopy(self)                
        if catCondition.startswith("N"):
            subDB.ndb(properN=properN, onlyProper=onlyProper)
        self._populateSubDB(catCondition, partition, subDB)
        subDB._ID += f'_{partition}-{catCondition}' 
        subDB._update.append(partition)
        if catCondition == "&":
            self._database += subDB._database
        elif partition == "pre":
            self._prenominalDomain = subDB
        elif partition == "post":
            self._postnominalDomain = subDB
        self._update.append((partition, catCondition))
            

    def _populateSubDB(self, catCondition, partition, subDB):
        """
        aux_@_partitionize


        Parameters
        ----------
        catCondition : TYPE
            DESCRIPTION.
        partition : TYPE
            DESCRIPTION.
        subDB : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        eliminate = []
        for np in subDB.database:
            if np.hasCat(catCondition):
                index = np.getIndex(catCondition)
            else:
                eliminate.append(np)
                continue
            if partition == "pre" and index > 0:
                startIdx = 0
                endIdx = index
                length = index
            elif partition == "post" and np.length - index > 1:
                startIdx = index + 1
                endIdx = np.length
                length = np.length - index - 1
            else:
                eliminate.append(np)
                continue
            np._ID = np.ID + f'_{partition}-{catCondition}'
            np._length = length
            np._lemmata = np._lemmata[startIdx:endIdx]
            np._pattern = {0: np.pattern[0][startIdx:endIdx],
                           1: np.pattern[1][startIdx:endIdx],
                           2: np.pattern[2][startIdx:endIdx],
                           3: np.pattern[3][startIdx:endIdx],
                           7: np.pattern[7][startIdx:endIdx]}
            np._categories = {0: set(cat for cat in np._pattern[0]),
                              1: set(cat for cat in np._pattern[1]),
                              2: set(cat for cat in np._pattern[2]),
                              3: set(cat for cat in np._pattern[3]),
                              7: set(cat for cat in np._pattern[7])}
        _ = [subDB._database.remove(np)
             for np in eliminate]
            

    def getPartitions(self, level=2, partition="pre", countFrom="left"):
        """
        returns for every category cat in the current DB the number (num) of 
        occurrences of 
        
        all categories c_2 
        
        returns NP-internal category distribution of the current DB ordered by 
        slots/columns; NPs are construed as aligned -- order="precede" <==> 
        left-aligned, order="follow" <==> right-aligned -- 
        
        "precede:"
        ||column_1  |column_2  |column_3  |column_4  | . . . . .|column_n  || 
        =====================================================================
        ||cat1: num |cat1: num |cat1: num |cat1: num |cat1: num |cat1: num ||
        ||cat2: num |cat2: num |cat2: num |cat2: num |cat2: num |......... ||
        ||cat3: num |cat3: num |cat3: num |cat3: num |......... |......... ||
        ||cat4: num |cat4: num |cat4: num |......... |......... |......... ||
        ||cat5: num |cat5: num |......... |......... |......... |......... ||
        ||cat6: num |......... |......... |......... |......... |......... ||
        
        "follow:"
        ||column_-1 |column_-2 |column_-3 |column_-4 | . . . . .|column_-n || 
        =====================================================================
        ||cat1: num |cat1: num |cat1: num |cat1: num |cat1: num |cat1: num ||
        ||......... |cat2: num |cat2: num |cat2: num |cat2: num |cat2: num ||
        ||......... |......... |cat3: num |cat3: num |cat3: num |cat3: num || 
        ||......... |......... |......... |cat4: num |cat4: num |cat4: num ||
        ||......... |......... |......... |......... |cat5: num |cat5: num ||
        ||......... |......... |......... |......... |......... |cat6: num ||
        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.
        partition : TYPE, optional
            DESCRIPTION. The default is "pre".
        countFrom : TYPE, optional
            DESCRIPTION. The default is "left".

        Raises
        ------
        ValueError
            DESCRIPTION:    if {direction} != "left" || "right"
        DataBaseNotInitializedError
            DESCRIPTION:    if partition has not been initialized, 
                            @see:  partitionize

        Returns
        -------
        TYPE:               tuple [ {Counter}, {Counter} ... ]
            DESCRIPTION:    each Counter object represents a column in the current partition; 
                            the count indicates the number of occurrences of categories 
                            in that column

        """
        if not countFrom.lower() in {"left", "right"}:
            raise ValueError("Unknown direction. Choose 'left' or 'right' !")        
        workingDB = self._getWorkingDB(partition)        
        maxCats = max(np.length 
                      for np in workingDB)
        out = [[] for i in range(maxCats)]        
        for np in workingDB:
            if countFrom.lower() == "left":
                startIdx, endIdx, steps = 0, np.length, 1
            else: 
                startIdx, endIdx, steps = -1, -np.length - 1, -1
    
            for idx in range(startIdx, endIdx, steps):
                out[idx].append(np.getPatt(level)[idx])
        out = [Counter(cats)
               for cats in out]
        return tuple(out)


    def _getWorkingDB(self, partition):
        if partition == "full":
            workingDB = self.database
        elif self.prenominalDomain.ID != "Old Empty-ese" and partition == "pre":
            workingDB = self.prenominalDomain.database
            
        elif self.postnominalDomain.ID != "Old Empty-ese" and partition == "post":
            workingDB = self.postnominalDomain.database
        else:
            msg = f'\nInitialize {partition}-partition <==>' \
                + f' partitionize("{partition}") !'
            raise DataBaseNotInitializedError(msg)
        return workingDB




###########################################################
###     Schrödinger's CATs 



    def _getCommonMultiple(self, numberOfColums):
        """
        aux_@_probabilize

        Parameters
        ----------
        numberOfColums : TYPE
            DESCRIPTION.

        Returns
        -------
        commonMultiple : TYPE
            DESCRIPTION.

        """
        commonMultiple = 1
        for i in range(2, numberOfColums+1):
            commonMultiple *=i
        return commonMultiple
        
    
    def _assignScores(self, patternLength, commonMultiple):
        """
        aux_@_probabilize

        Parameters
        ----------
        patternLength : TYPE
            DESCRIPTION.
        commonMultiple : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        scores = []
        step = commonMultiple // (patternLength - 1)
        score = commonMultiple
        for slots in range(patternLength):
            scores.append(score)
            score -= step
        scores[0] -= 0.1
        scores[-1] += 0.1
        return tuple(scores)
        
    def _evaluatePattern(self, pattern, commonMultiple):
        """
        aux_@_probabilize        

        Parameters
        ----------
        pattern : TYPE
            DESCRIPTION.
        commonMultiple : TYPE
            DESCRIPTION.

        Returns
        -------
        dict
            DESCRIPTION.

        """
        scores = self._assignScores(len(pattern), commonMultiple)
        return {pattern[i]: scores[i]
                for i in range(len(pattern))}
        

    def _getScores(self, level, commonMultiple, score_toCat, catCount):
        """
        aux_@_probabilize        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.
        commonMultiple : TYPE
            DESCRIPTION.
        score_toCat : TYPE
            DESCRIPTION.
        catCount : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for patt in self.patternize(level):
            if len(patt) < 2:
                continue
            catCount += Counter(patt)
            scoresAsigned = self._evaluatePattern(patt, commonMultiple)
            for cat in scoresAsigned:
                score_toCat[cat] += scoresAsigned[cat]


    def probabilize(self, level=2):
        """
        Calculates "Distance from N" (negative numbers in the prenominal domain);
        possible maximum = number of max categories / columns in the respective domain.

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        _y_factor = 1
        if "pre" in self.ID: 
            _y_factor = (-1)
        catCount = Counter()
        commonMultiple = self._getCommonMultiple(self.maxPatternSize)
        re_factor = commonMultiple // self.maxPatternSize
        score_ofCat = {cat : 0 
                       for cat in self.categorize(level)}
        self._getScores(level, commonMultiple, score_ofCat, catCount)
        eliminate = [cat 
                     for cat in score_ofCat
                     if not score_ofCat[cat]]
        _ = [score_ofCat.pop(cat)
             for cat in eliminate]
        for cat in score_ofCat:
            score_ofCat[cat] = score_ofCat[cat] / catCount[cat]
            finalScore = score_ofCat[cat] / re_factor
            finalScore *= _y_factor
            if _y_factor == 1: 
                finalScore = self.maxPatternSize - finalScore
            finalScore = round(finalScore, 2)
            if not finalScore and "-N" in self.ID:      
                finalScore = 0.01 * _y_factor
            score_ofCat[cat] = finalScore
        return Counter(score_ofCat)
                 
            
        
        

        
        
        
    def _cats_ge_pair(self, level=2):
        """
        auxiliary method

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return Counter([cat
                        for patt in self._patterns_ge_pair(level=level)
                        for cat in patt])
        
        
        
    def _patterns_ge_pair(self, level=2):
        """
        auxiliary method        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        list
            DESCRIPTION.

        """
        patterns = self.patternize(level=level)
        return [list(pattern) 
                for pattern in patterns
                for times in range(patterns[pattern])
                if len(pattern) > 1]
            


            

    def pair_wize(self, level=2, distance=1, order="precede"):
        """
        Calculates cooccurrences for any two categories at {level}; 
        {distance} indicates the max distance between two categories; z.B.
        distance=1 (immediate adjacency): 
            - cat1, cat2    --> (cat1, cat2) is a pair,
            - cat1, X, cat2 --> (cat1, cat2) is not a pair bc distance would be 2
        distance=2:
            - cat1, cat2    --> (cat1, cat2) is a pair,
            - cat1, X, cat2 --> (cat1, cat2) is a pair  
        Returns a dictionary of type Cat: Counter( {cat1, cat2 ...} ) 
        where the relation between Cat and cat1, cat2 .. catn is specified 
        by the parameter {order}, i.e. "precede" or "follow". 
        
        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.
        distance : TYPE, optional
            DESCRIPTION. The default is 1.
        order :     TYPE, optional
                    possible values: 
                        - "precede",
                        - "follow"
                    The default is "precede".

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        out : TYPE
            DESCRIPTION.

        """
        if distance < 1 or distance >= self.maxPatternSize:
            msg = f'Distance value {distance} is inapproriate!\n' \
                + f'Possible values x: 1 <= x <= {self.maxPatternSize - 1}'
            raise ValueError(msg)
        patterns = self._patterns_ge_pair(level=level)
        cats = self._cats_ge_pair(level=level)
        
        out = dict()
        for firstCat in cats:
            otherCats = []
            for secondCat in cats: 
                self._addIf_PairMatch(firstCat, secondCat, patterns, distance, 
                                  order, otherCats)
            if otherCats:
                out[firstCat] = Counter(otherCats)
        return out
                    



    def _addIf_PairMatch(self, firstCat, secondCat, patterns, distance, order, 
                     otherCats):
        """
        auxiliary method        


        Parameters
        ----------
        firstCat : TYPE
            DESCRIPTION.
        secondCat : TYPE
            DESCRIPTION.
        patterns : TYPE
            DESCRIPTION.
        distance : TYPE
            DESCRIPTION.
        order : TYPE
            DESCRIPTION.
        otherCats : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for patt in patterns:
            pattern = deepcopy(patt)
            firstMatch = False
            if order == "follow":
                pattern.reverse()
            if firstCat in pattern:
                preIdx = pattern.index(firstCat)
                firstMatch = True
            if firstMatch and len(pattern) - preIdx > distance and \
                pattern[preIdx+distance] == secondCat:
                    otherCats.append(secondCat) 




                      #  specfics for 2 cats
    def conCatenize(self, cat1, cat2, level=2):
        """
        Tests / measures the "pair-hood" relations for two given categories, 
        i.e. how often 
            - do {cat1} and {cat2} cooccur 
            - does {cat1} precede {cat2}, 
            - does {cat2} precede {cat1}?
        Returns a dictionary.
        

        Parameters
        ----------
        cat1 : TYPE
            DESCRIPTION.
        cat2 : TYPE
            DESCRIPTION.
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        dict
            DESCRIPTION.

        """
        
        
        patterns = self._patterns_ge_pair(level=level)
        cooccurrences = 0
        cat1_precedes = 0
        cat2_precedes = 0
        cat1_score = 0
        cat2_score = 0
        for pattern in patterns:
            if cat1 in pattern and cat2 in pattern:
                cooccurrences += 1
                Idx1 = pattern.index(cat1)
                Idx2 = pattern.index(cat2)
                if Idx1 < Idx2:
                    cat1_precedes += 1
                else:
                    cat2_precedes += 1
                cat1_score += (Idx2 - Idx1)
                cat2_score += (Idx1 - Idx2)
                if pattern.count(cat1) > 1:
                    cooccurrences += 1
                    Idx1 = pattern.index(cat1, Idx1+1)
                    if Idx1 < Idx2:
                        cat1_precedes += 1
                    else:
                        cat2_precedes += 1
                        cat1_score += (Idx2 - Idx1)
                        cat2_score += (Idx1 - Idx2)
                if pattern.count(cat2) > 1:
                    cooccurrences += 1
                    Idx2 = pattern.index(cat2, Idx2+1)
                    if Idx1 < Idx2:
                        cat1_precedes += 1
                    else:
                        cat2_precedes += 1
                        cat1_score += (Idx2 - Idx1)
                        cat2_score += (Idx1 - Idx2)
        return {'cooccur': cooccurrences,
                f'{cat1}' : {"precedes" : cat1_precedes,
                             "score" : cat1_score}, 
                f'{cat2}' : {"precedes" : cat2_precedes,
                             "score" : cat2_score} 
                }
                    


    ########### TO DO !!! ###############
    def _replace_CatLabelsList(self, labels):
        for i in range(len(labels)):
            if labels[i] in self._replaceCat:
                labels[i] = self._replaceCat[labels[i]]
        print(labels)
        return labels
                
    ########### TO DO !!! ###############
    def _replace_CatLabelsInd(self, cat):
        if cat in self._replaceCat:
            return self._replaceCat[cat]
        else:
            return cat
                
    ########### TO DO !!! ###############
    


    def showPobabilisticCatDistribution(self, level=2):
        """
        NB: TO DO -- this method is still under construction and not fully functional;
                      but can be called (--> preliminary view).

        Parameters
        ----------
        level : TYPE, optional
            DESCRIPTION. The default is 2.

        Returns
        -------
        None.

        """
        import matplotlib.pyplot as plt
        import networkx as nx
        diGraph = nx.DiGraph()
        _cats = {cat : self._cats_ge_pair(level=level)[cat]
                 for cat in self._cats_ge_pair(level=level)
                 if self._cats_ge_pair(level=level)[cat] > 10 }
        node_labels, positions = self._get_labelledCatPositions(level, _cats)            
        diGraph.add_nodes_from(node_labels)
        edges, edge_colours = self._getEdges(level, _cats)
        diGraph.add_edges_from(edges)

        colour_labels = self._catColours(node_labels)
        
        
        nx.draw(diGraph,
                with_labels=True,
                
                node_color=colour_labels,
                node_size=1000, 
                pos=positions,
                edge_color=edge_colours,
                width=1,       # width of edges
                
                font_color="black",
                font_size=10,
                font_family="Times New Roman")

        plt.margins(0.2)
        plt.grid()
        plt.show()
            

    def _replace_CatLabelsDict(self, labels):
        out = dict()
        for cat in labels:
            if cat in self._replaceCat:
                out[self._replaceCat[cat]] = labels[cat]
            else: 
                out[cat] = labels[cat]
        return out
                
                

    def _get_labelledCatPositions(self, level, _cats):
        """
        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.
        _cats : TYPE
            DESCRIPTION.

        Returns
        -------
        labels : TYPE
            DESCRIPTION.
        positions : TYPE
            DESCRIPTION.

        """
        _catProbabilities = self.probabilize(level=level)
        cat_probScores = {cat : _catProbabilities[cat] 
                          for cat in _cats} 
        catCount = {cat : _cats[cat] 
                    for cat in _cats}
        labels = []
        positions = dict()
        for cat in _cats:
            
            if cat == "Md":         ##
                continue
            
            labels.append(cat)
            positions[cat] = (cat_probScores[cat], catCount[cat])
            
      #      positions = self._replace_CatLabelsDict(positions)
      #      labels = self._replace_CatLabelsList(labels)

        return labels, positions 
            
            
    def _catColours(self, labels):
        colours = []
        for cat in labels:
            if cat in {"Mod", "Adj", "ALex", "AFun", "GenP", "RC", "CC", "Mdmd",
                       "Md.Aj", "Md.Aj.Lx", "Md.Aj.Fn", "Md.Nu/WQ", 
                       "Md.Nu/WQ.Nu","Md.Nu/WQ.WQ"}:
                colours.append("fuchsia")
            elif cat in {"Dem", "H", "Poss"}:
                colours.append("springgreen")
            elif cat in {"Q", "WQ", "Num"}:
                colours.append("cyan")
            else:
                colours.append("yellow")
        return colours
                
            
        
        
    def _getEdges(self, level, _cats):
        """
        

        Parameters
        ----------
        level : TYPE:     int
                          The default is 2.
            DESCRIPTION:  level of sub-categorization.
        _cats : TYPE
            DESCRIPTION.

        Returns
        -------
        edges : TYPE
            DESCRIPTION.
        edge_colours : TYPE
            DESCRIPTION.

        """
        catPairs = combinations(_cats, 2)
        edges = []        
        edge_colours = []        
        try:
            while True:
                pair = next(catPairs)
                cat1, cat2 = pair[0], pair[1]
                edgeTest = self.conCatenize(cat1, cat2, level=level)
                
                if cat1 == "Md" or cat2 == "Md":    ###
                    continue
                
                if edgeTest["cooccur"] < 20:    #choose a minimum
                    continue
                if edgeTest["cooccur"] > 100:
                    edge_colours.append("red")
                elif 100 >= edgeTest["cooccur"] > 55:
                    edge_colours.append("green")
                elif 55 >= edgeTest["cooccur"] >= 20:
                    edge_colours.append("blue")
                else: 
                    edge_colours.append("gray")

                if edgeTest[cat1]["precedes"] > edgeTest[cat2]["precedes"]:
                    preCat, postCat = cat1, cat2
                else: 
                    preCat, postCat = cat2, cat1
                
                edges.append((preCat, postCat)) 
        except StopIteration:
            return edges, edge_colours

        
        
        
        
                    

###################################################################        


############################################
###     DUNDER


    def __repr__(self):
      return "<class \'Patternize\'>"  
        
    def __str__(self):
        return f'[Patternize::Database; Language: {self.ID}; Size: {self.size}]'

    def __len__(self):
        return self.size








###############################################################################
###############################################################################
###
###     PATTERN (class)
###

class Pattern:
    
    def __init__(self, np):
        """
        ==> for internal use only!

        Creates a Pattern object <==> NP construed as .database element 
                                        @ Patternize:Database
        

        Parameters
        ----------
        np : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self._ID = np["DB_item_id"]
        self._length = np["Length"]
        self._lemmata = np["Lemmata"]
        self._pattern = {0: tuple(np["Patterns"]["cat0"]),
                         1: tuple(np["Patterns"]["cat1"]),
                         2: tuple(np["Patterns"]["cat2"]),
                         3: tuple(np["Patterns"]["cat3"])}
        self._categories = {0: set(np["Categories"]["cat0"]),
                            1: set(np["Categories"]["cat1"]),
                            2: set(np["Categories"]["cat2"]),
                            3: set(np["Categories"]["cat3"])}

        self._macroPattern = tuple()
        self._allCats = set()
        

    
    @property
    def ID(self):
        """
        -- NPEGL ID of this NP indicating the language and id nr.

        """
        return self._ID
    
    @property
    def length(self):
        """
        -- Length of this NP measured in numbers of categories

        """
        return self._length
    
    @property 
    def lemmata(self):
        """
        -- Lexical items involved in this NP

        """
        return tuple(self._lemmata)
    
    @property
    def pattern(self):
        """
        -- Patterns instantiated by this NP; 
            is returned as dictionary and can be accessed via >> .pattern[level] <<

        """
        return self._pattern
        
    def getPatt(self, level):
        """
        Returns the pattern at {level} instantiated by this NP
        Parameters
        ----------
        level : TYPE:     int
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        TYPE:               tuple of str        
            DESCRIPTION:    sequence of categories <==> pattern

        """
        return self.pattern[level]
        
    @property
    def categories(self):
        return self._categories
    
    @property
    def allCats(self):
        self._allCats = self._categories[0].union(self._categories[1])
        self._allCats = self._allCats.union(self._categories[2])
        self._allCats = self._allCats.union(self._categories[3])
        return self._allCats
        

    @property
    def macroPattern(self):
        self._macroPattern = []
        for position in range(self.length):
            cats = {self._pattern[catLevel][position]
                    for catLevel in self._pattern}
            cats.add(self.lemmata[position])
            self._macroPattern.append(cats)
        return tuple(self._macroPattern)

        
    def getCat(self, level):
        """
        Parameters
        ----------
        level : TYPE:     int
            DESCRIPTION:  level of sub-categorization.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.categories[level]  
    

    def hasCat(self, cat):
        """
        Diagnostic method: does self have {cat}?

        Parameters
        ----------
        cat : TYPE:       str
            DESCRIPTION:  (sub-)category label

        Returns
        -------
            TYPE:         bool
            DESCRIPTION:  True if self contains {cat}

        """
        return cat in self.allCats

    
    def getIndex(self, cat):
        """
        Returns the index position (of the first occurrence) of {cat} 
        within the pattern -- iff present, else -1 if {cat} is not present in the pattern.

        Parameters
        ----------
        cat : TYPE
            DESCRIPTION.

        Returns
        -------
        index : TYPE
            DESCRIPTION.

        """
        index = -1
        for cats in self.macroPattern:
            if cat in cats:
                index = self.macroPattern.index(cats)
        return index
    
    
        
            
    
    def _replace_CatLabel(self, cat):
        if cat in self._replaceCat:
            return self._replaceCat[cat]
        else:
            return cat
    

    def __repr__(self):
      return "<class \'Pattern\'>"  
        
    def __str__(self):
        return f'[Patternize::Pattern; @ID: {self.ID}; \n patt2: {self.getPatt(2)}]'



       
###############################################################################
###############################################################################
###
###     COMBFLEX (auxiliary class)
###
        
class CombFlex:
    
    def __init__(self, combination, length, sm_pattern, align, 
                 pattern_threshold, group_threshold, count, groupCount,
                 permutations, catCondition):
        """
        ==> for internal use only!
        
        Implementation of 'Combinatorial Flexibility', 
        -- this data object is returned by Patternize.combinatorialFlexibility

        Parameters
        ----------
        combination : TYPE
            DESCRIPTION.
        length : TYPE
            DESCRIPTION.
        sm_pattern : TYPE
            DESCRIPTION.
        align : TYPE
            DESCRIPTION.
        pattern_threshold : TYPE
            DESCRIPTION.
        group_threshold : TYPE
            DESCRIPTION.
        count : TYPE
            DESCRIPTION.
        groupCount : TYPE
            DESCRIPTION.
        permutations : TYPE
            DESCRIPTION.
        catCondition : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        self._combination = combination
        self._length = length
        self._sm_pattern = sm_pattern
        self._alignment = align
        self._pattern_threshold = pattern_threshold
        self._group_threshold = group_threshold
        self._count = count
        self._groupCount = groupCount
        self._catCondition = catCondition
        self._permutations = permutations
        
    @property
    def combination(self):
        return set(self._combination)
    
    @property
    def length(self):
        return self._length
    
    @property
    def sm_pattern(self):
        return self._sm_pattern
        
    @property
    def alignment(self):
        return str(self._alignment)
    
    @property
    def pattern_threshold(self):
        return self._pattern_threshold

    @property
    def group_threshold(self):
        return self._group_threshold

    @property
    def count(self):
        return self._count
    
    @property
    def groupCount(self):
        return self._groupCount
    
    @property
    def catCondition(self):
        return self._catCondition


    @property
    def permutations(self):
        return self._permutations
    
    
    def getAttestations(self):
        return {per_mutation : self.permutations[per_mutation]["Attestations"]
                for per_mutation in self.permutations.keys()}
        
        
    def showAttestations(self):
        for perm in self.getAttestations():
            print(perm)
            for patts in self.getAttestations()[perm]:
                print("          ", patts.getPatt(2))
            print(" = " * 17)
            print()
    
    def info(self):
        return f'[Patternize::CombFlex; @combination: {self.combination};\n' \
            + f' sm_pattern: {self._sm_pattern}; alignment: {self.alignment}]'

    def getCombFlex(self):
        att =  [self.permutations[p]['Count'] 
                for p in self.permutations.keys()
                if self.permutations[p]['Count'] != False] 
        out = str(len(att)) + ' / ' +  str(len(self.permutations))
        return out


    def mirror(self):
        if self._length != 3:
            msg = 'Pattern must be of length 3!'
            raise IllegalValueException(msg)
        if (self._catCondition != "N"  and self._catCondition != "N.C" and 
            self._catCondition != "N.P"):
            msg = 'A nominal category must be present!'
            raise NoNounException(msg)
        cats = list(self._combination)
        cats.remove(self._catCondition)
        cat1 = cats[0]
        cat2 = cats[1]
        outStr = "".join(["\n  Special Ordering -- Mirror order:\n",
              f"\t {cat1} >> {cat2}  {self._catCondition}: ",
              f"{self._permutations[(cat1, cat2, self._catCondition)]['Count']:>4}",
              f"\t\t\t{self._catCondition}  {cat2} << {cat1}: ",
              f"{self._permutations[(self._catCondition, cat2, cat1)]['Count']:>4}\n",
              f"\t {cat2} >> {cat1}  {self._catCondition}: ",
              f"{self._permutations[(cat2, cat1, self._catCondition)]['Count']:>4}",
              f"\t\t\t{self._catCondition}  {cat1} << {cat2}: ",
              f"{self._permutations[(self._catCondition, cat1, cat2)]['Count']:>4} "])
        toStr = str(self)[:-67]
        return toStr + outStr + "\n\n = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n"


    def precedence(self):
        if self._length != 3:
            msg = 'Pattern must be of length 3!'
            raise IllegalValueException(msg)
        if (self._catCondition != "N"  and self._catCondition != "N.C" and 
            self._catCondition != "N.P"):
            msg = 'A nominal category must be present!'
            raise NoNounException(msg)
        cats = list(self._combination)
        cats.remove(self._catCondition)
        cat1 = cats[0]
        cat2 = cats[1]
        outStr = "".join(["\n  Special Ordering -- Precedence:\n",
              f"\t {cat1} >> {cat2}  {self._catCondition}: ",
              f"{self._permutations[(cat1, cat2, self._catCondition)]['Count']:>4}",
              f"\t\t\t{self._catCondition}  {cat1} >> {cat2}: ",
              f"{self._permutations[(self._catCondition, cat1, cat2)]['Count']:>4}\n",
              f"\t {cat2} >> {cat1}  {self._catCondition}: ",
              f"{self._permutations[(cat2, cat1, self._catCondition)]['Count']:>4}",
              f"\t\t\t{self._catCondition}  {cat2} >> {cat1}: ",
              f"{self._permutations[(self._catCondition, cat2, cat1)]['Count']:>4} "])
        toStr = str(self)[:-67]
        return toStr + outStr + "\n\n = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n"



    def flanked(self):
        if self._length != 3:
            msg = 'Pattern must be of length 3!'
            raise IllegalValueException(msg)
        if (self._catCondition != "N"  and self._catCondition != "N.C" and 
            self._catCondition != "N.P"):
            msg = 'A nominal category must be present!'
            raise NoNounException(msg)
        cats = list(self._combination)
        cats.remove(self._catCondition)
        cat1 = cats[0]
        cat2 = cats[1]
        outStr = "".join(["\n  Special Ordering -- Flanked:\n",
              f"\t {cat1}  >>  {self._catCondition}  >>  {cat2}: ",
              f"{self._permutations[(cat1, self._catCondition, cat2)]['Count']:>4}",
              f"\t\t\t{cat2}  >>  {self._catCondition}  >>  {cat1}: ",
              f"{self._permutations[(cat2, self._catCondition, cat1)]['Count']:>4}\n"])
        toStr = str(self)[:-67]
        return toStr + outStr + "\n\n = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n"







    def __str__(self):
        permCount = [(perm, self.permutations[perm]["Count"]) 
                     for perm in self.permutations ]         
        outStr = f'\nCombination: {self.combination}\n\n'                   \
            + f'  SM_Pattern: {self.sm_pattern}\n\n'                        \
            + f'  Combinatorial Flexibility: {self.getCombFlex()}  \n\n'                          \
            + f'  Alignment: {self.alignment}\n\n'                          \
            + f'  GroupCount: {self.groupCount}\n\n'                        \
            + f'  Permutations:   {permCount[0][0]}: {permCount[0][1]} \n'  \
                
        for perm in permCount[1:]:
            outStr += f'                  {perm[0]}: {perm[1]}\n'
        outStr = outStr +  "\n\n = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n"           
        return outStr
        
    def __repr__(self):
        return f"<class \'CombFlex; Combination: {self.combination}\'>" 







       









       
###############################################################################
###############################################################################
###
###     EXCEPTIONS
###


class DataBaseNotInitializedError(Exception):
    def __init__(self, message):
        super().__init__(message)


class InvalidOrderException(Exception):
    def __init__(self, message):
        super().__init__(message)



class IllegalValueException(Exception):
    def __init__(self, message):
        super().__init__(message)


class NoNounException(Exception):
    def __init__(self, message):
        super().__init__(message)













def main():
    
    print("\n***********************************************")
    print("********* Patternization@NPEGL_v.9 ************")
    print("***********************************************\n")


    print("===============================================")
    print("Language key (default <=> Old Icelandic):\n") 
    print("1. Old Icelandic: 'oice' ")
    print("2. Old Saxon: 'osax' ")
    print("3. Old English: 'oeng' ")
    print("===============================================\n")


 

main()


cats = ['Md.Nu/WQ.WQ', 'Md.Aj.Fn', 'Md.Aj.Fn', 'Md.Aj.Lx', 'Md.Aj.Lx', 
        'Poss', 'Dem', 'Q', "N.C"]

catSelection = [
        'Md.Nu/WQ.WQ', 'Md.Aj.Fn', 'Md.Aj.Fn', 'Md.Pos', 'Dgcm.Br', "N.C",
        'Md.Nu/WQ.Nu', 'GenP', 'Q', 'Dem', 'Mdmd', 'Mdcm.N', 'CC.Fi',
        'H', 'Mdcm', 'RC', 'Md.Nu/WQ', 'App', 'Per', 'IXP', 'Dgcm', 
        'Poss', 'Dgcm.Mk', 'Md.Aj', 'Md.Aj', 'disco', 'CC', 'N.C', 'PP',
        'Md.Aj.Lx', 'Md.Aj.Lx', 'Mdcm.P', 'Assoc', 'CC.Nf', 'Md', 'Adv']



