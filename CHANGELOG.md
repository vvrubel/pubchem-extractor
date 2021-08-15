# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
* changelog
### Changed
* project structure: source code includes three modules – api, cli and db
* project dependencies
* mypy config

## [0.1.0] - 2021-08-08
### Added
* a command-line interface for fetching Pubchem data and uploading it to your local MongoDB server;
* API to interact with local MongoDB data where you can run a substructure search and get the results in two different ways: by page or in the summary form.
