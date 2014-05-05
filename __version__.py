"""This identifies the version of heatsource. We use the format specified in PEP 386"""

# major 
#   0 - 99999
#   Major designates a major revision. This occurs when you are adding a lot of features,
#   breaking backward-compatibility, or drastically changing the model design.

# minor
#   0 - 99999
#   Minor usually groups moderate changes to the code like bug fixes and feature improvements that 
#   may impact model operation. End users can upgrade with no backward compatibility issues.
#   In some cases where major new features are added or others are planned to be removed, end users will be
#   notified with deprecation warnings but backward-compatibility is maintained.

# micro
#   0 - 99999
#   Micro changes are devoted to minor bug fixes that have minimal impact on model operation

# prerel (Pre Releases)
#   'a' = alpha, 'b' = beta, 'c' = release candidate, '' = final

#   'a' = alpha, 
#   Early pre-releases. A lot of changes can occur between alphas and the
#   final release, like feature additions or refactorings. But they are minor changes and the
#   software should stay pretty unchanged by the time the first beta is reached.

#   'b' = beta
#   No new features are added and developers are tracking remaining bugs.

#   'c' = release candidate
#   A release candidate is an ultimate release before the final release. 
#   Unless something bad happens, nothing is changed.

# prerelversion
#    1 - 99999, '' = final
#    pre release version number

major = 9
minor = 0
micro = 0
prerel = 'b'
prerelversion = '7' 
version_info = (major,minor,micro,prerel,prerelversion)
version_string = "%s.%s.%s%s%s" % version_info
