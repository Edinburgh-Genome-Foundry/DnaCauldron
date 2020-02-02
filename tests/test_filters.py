import dnacauldron as dc

def test_filters():
    to_record = dc.sequence_to_biopython_record

    filter1 = dc.NoRestrictionSiteFilter("BsmBI")
    assert filter1(to_record("ATGATGATG"))
    assert not filter1(to_record("ACGTCTCTTG"))
    assert not filter1(to_record("ACGGAGACGG"))

    filter2 = dc.TextSearchFilter("GFP", is_forbidden=True)
    record = to_record("ATCGCGTGCGTGCACCACACGT")
    assert filter2(record)
    dc.annotate_record(record, location=(20, 40), label="Here's some GFP!")
    assert not filter2(record)