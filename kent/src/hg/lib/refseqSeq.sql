# used in KG build.
CREATE TABLE refseqSeq (
  name varchar(255) NOT NULL default '',
  seq longblob NOT NULL,
  PRIMARY KEY  (name(32))
) TYPE=MyISAM;

