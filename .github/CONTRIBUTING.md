# Contribution Guidelines for CMS-HGCAL

## Pull-requests

When creating a pull-request, state clearly what it is about by choosing a
meaningful title.

Give a bulleted list of things you address (mentioning related issues).

We are happy if you create pull-requests also if you feature is not ready yet.
Please mark them as such by adding `[WIP]` to the beginning of the title. The purpose
of this is, for example, that you want to let other people know you are working
on a given issue or to get early feedback on your code. If you would like to receive feedback,
mention it in the description. By default we do not review pull requests marked
as WIP. We propose to have a check list of things that still need to be done
using the markdown format `- [ ]` (especially if you want feedback).

## Issue tracking

Use the [github issue tracker](https://github.com/CMS-HGCAL/reco-ntuples/issues). Reference
the issues that you are working on. If you notice a problem / bug, consider first
creating an issue and then referring to it in your pull-request and commit
messages with `#[issue-id]`.

## Code style

Please do not change indentation, but follow the indentation currently used.
Consider using a `clang-format` plugin for your favourite editor or run `clang-format` after your changes to fix the format automatically.
