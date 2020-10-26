<template>
  <div class="notification is-success is-light" v-show="!isHidden">
    <button class="delete" @click="close"></button>
    {{ banner_text }}
    Learn more <a :href="banner_learn_more_link">here</a>
  </div>
</template>

<script>
export default {
  name: "CookieBanner",
  computed: {
    isHidden() {
      return this.getState();
    },
  },
  data() {
    return { isHiddenState: false };
  },
  methods: {
    close() {
      // to prevent the node build from failing since
      // localStorage is not available in node
      if (this.remember && typeof localStorage !== `undefined`) {
        localStorage.setItem("cookie-banner", JSON.stringify(true));
      }
      return (this.isHiddenState = true);
    },
    getState() {
      // to prevent the node build from failing since
      // localStorage is not available in node
      const savedState =
        typeof localStorage !== `undefined`
          ? localStorage.getItem("cookie-banner")
          : this.remember;
      return savedState && this.remember
        ? JSON.parse(savedState)
        : this.isHiddenState;
    },
  },
  props: {
    banner_text: String,
    banner_learn_more_link: String,
    // by default, cookie banner will remain closed
    // per user session if user have closed it once
    remember: {
      default: true,
      type: Boolean,
    },
  },
};
</script>

<style lang="sass" scoped>
@import "../../node_modules/bulma/bulma.sass"
</style>
