from genedesign.checkers.gc_content_checker import GCContentChecker

def test_gc_content_checker_valid():
    checker = GCContentChecker()
    # 50% GC content (4 G/C, 4 A/T) - This is biologically ideal!
    sequence = "ATGCATGC" 
    passes, gc_content = checker.run(sequence)
    
    assert passes == True
    assert gc_content == 0.5

def test_gc_content_checker_invalid_high():
    checker = GCContentChecker()
    # 100% GC content - Too high, should fail!
    sequence = "GCGCGCGC"
    passes, gc_content = checker.run(sequence)
    
    assert passes == False
    assert gc_content == 1.0
    
def test_gc_content_checker_invalid_low():
    checker = GCContentChecker()
    # 0% GC content - Too low, should fail!
    sequence = "ATATATAT"
    passes, gc_content = checker.run(sequence)
    
    assert passes == False
    assert gc_content == 0.0