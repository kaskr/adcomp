# Automatic deploy documentation

Setting up travis to deploy documentation on every push to master involves some steps that are easy to forget:

1. Generate ssh keys (private and public).
   https://help.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
2. Go to https://github.com/kaskr/adcomp/settings/keys and add **public** key. **Tick write access** .
3. Install travis-ci tools: `sudo apt install ruby ruby-dev` and `sudo gem install travis`.
4. Copy **private** key to `adcomp/bin/deploy_key` and run `travis encrypt-file deploy_key`. *Note the output* !
5. Edit `bin/deploy.sh` with the output from (4) i.e. 'openssl bla bla'.
