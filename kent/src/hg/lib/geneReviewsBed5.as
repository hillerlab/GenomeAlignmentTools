table geneReviewsBed5
"GeneReviews  bigBed 4 + with extra field for detail page"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Short Name of item"
    lstring description; "geneReviews item details in html"
    )
